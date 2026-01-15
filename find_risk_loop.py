"""
Identify GWAS risk-variant–anchored chromatin loops and link them to gene promoters (±2 kb of TSS).

Pipeline (high-level):
1) Load Mustache loops and compute loop distance for cis filtering.
2) Build loop anchors (A/B) as genomic intervals.
3) Intersect anchors with GWAS significant risk-increasing SNPs.
4) For each SNP-hit anchor, retrieve the opposite anchor (other side of the same loop).
5) Intersect opposite anchors with gene promoter windows (±2 kb around TSS) from GENCODE.
6) Output SNP → loop → gene links.

Outputs:
- risk_loop.csv
"""

import pandas as pd
import numpy as np
import pyranges as pr

# -----------------------------
# Parameters
# -----------------------------
MIN_DIST = 20_000        # Remove loops closer than 20 kb (cis only)
MAX_DIST = 5_000_000     # Upper bound for cis distance (adjust as needed)
GWAS_P_THRESH = 5e-8     # Genome-wide significance threshold
GWAS_BETA_SIGN = "pos"   # Keep beta > 0 (risk-increasing alleles)
PROMOTER_HALF_WINDOW = 2000  # ±2 kb around TSS

# -----------------------------
# 1) Load loops and filter by distance
# -----------------------------
loops = pd.read_table("mustache_SCALE_CD4.tsv")
loops.columns = [
    "chrA", "A_start", "A_end",
    "chrB", "B_start", "B_end",
    "FDR", "DETECTION_SCALE", "loop_id"
]

# Compute midpoint distance only for intra-chromosomal loops
midA = (loops["A_start"] + loops["A_end"]) // 2
midB = (loops["B_start"] + loops["B_end"]) // 2
same_chr = loops["chrA"] == loops["chrB"]
loops["distance"] = np.where(same_chr, (midA - midB).abs(), np.nan)

# Keep all trans loops; for cis loops apply distance bounds
loops_f = loops[
    (~same_chr) | ((loops["distance"] >= MIN_DIST) & (loops["distance"] <= MAX_DIST))
].copy()

# -----------------------------
# 2) Build anchors (A and B) as intervals
# -----------------------------
A = loops_f[["chrA", "A_start", "A_end", "loop_id"]].copy()
A.columns = ["chr", "start", "end", "loop_id"]
A["side"] = "A"

B = loops_f[["chrB", "B_start", "B_end", "loop_id"]].copy()
B.columns = ["chr", "start", "end", "loop_id"]
B["side"] = "B"

anchors = pd.concat([A, B], ignore_index=True)

# Convert anchors to PyRanges (requires Chromosome/Start/End)
anchors_gr = pr.PyRanges(
    anchors.rename(columns={"chr": "Chromosome", "start": "Start", "end": "End"})
           [["Chromosome", "Start", "End", "loop_id", "side"]]
)

# Build an index to quickly fetch the opposite anchor for a given (loop_id, side)
anchor_idx = anchors.set_index(["loop_id", "side"])

# -----------------------------
# 3) Load GWAS SNPs and keep significant risk-increasing variants
# -----------------------------
ma = pd.read_csv("./GWAS_summary/NG_2021_JT.ma")

# Filter by significance and effect direction (beta > 0)
ma = ma[(ma["p_value"] < GWAS_P_THRESH) & (ma["beta"] > 0)]

# Build SNP intervals: [pos, pos+1)
snps_df = ma[["chromosome", "base_pair_location", "rsID"]].copy()
snps_df.columns = ["Chromosome", "Start", "rsID"]
snps_df["End"] = snps_df["Start"] + 1
snps_gr = pr.PyRanges(snps_df[["Chromosome", "Start", "End", "rsID"]])

# -----------------------------
# 4) Intersect anchors with GWAS SNPs (anchor hits)
# -----------------------------
hits_gr = anchors_gr.join(snps_gr)  # interval overlap join
snp_hits = hits_gr.as_df()

# Standardize column names from PyRanges output
# Left (anchors): Chromosome, Start, End
# Right (SNPs): Chromosome_b, Start_b, End_b, rsID
snp_hits = snp_hits.rename(columns={
    "Chromosome": "chr",
    "Start": "start",
    "End": "end",
    "Start_b": "snp_start",
    "End_b": "snp_end"
})[["chr", "start", "end", "loop_id", "side", "snp_start", "snp_end", "rsID"]]

# -----------------------------
# 5) For each SNP-hit anchor, retrieve the opposite anchor in the same loop
# -----------------------------
opp_side = {"A": "B", "B": "A"}

def fetch_opposite_anchor(row: pd.Series) -> pd.Series:
    """
    Given a row with (loop_id, side), fetch the opposite anchor coordinates.
    """
    oside = opp_side[row["side"]]
    d = anchor_idx.loc[(row["loop_id"], oside)]
    return pd.Series(
        [d["chr"], d["start"], d["end"], oside],
        index=["opp_chr", "opp_start", "opp_end", "opp_side"]
    )

snp_hits = pd.concat([snp_hits, snp_hits.apply(fetch_opposite_anchor, axis=1)], axis=1)

# Deduplicate opposite anchors before intersecting with promoters
opp_df = (
    snp_hits[["opp_chr", "opp_start", "opp_end", "loop_id"]]
    .dropna()
    .drop_duplicates()
    .rename(columns={"opp_chr": "Chromosome", "opp_start": "Start", "opp_end": "End"})
)

opp_gr = pr.PyRanges(opp_df[["Chromosome", "Start", "End", "loop_id"]])

# -----------------------------
# 6) Build gene promoter intervals (±2 kb around TSS) from GENCODE pickle
# -----------------------------
gtf = pd.read_pickle("gencode_gtf.pickle")

# Use protein-coding genes; compute TSS based on strand
genes = gtf[(gtf["feature"] == "gene") & (gtf["gene_type"] == "protein_coding")].copy()
genes["gene_id"] = genes["gene_id"].str.split(".").str[0]

# TSS: start for '+' strand, end for '-' strand
genes["TSS"] = np.where(genes["strand"] == "+", genes["start"].astype(int), genes["end"].astype(int))

# Promoter windows: ±2 kb around TSS
genes["Start"] = genes["TSS"] - PROMOTER_HALF_WINDOW
genes["End"] = genes["TSS"] + PROMOTER_HALF_WINDOW

# Keep autosomes chr1–chr22; normalize chromosome naming if needed
# Here the code expects chromosomes to be "1".."22" (no "chr" prefix) to match your loops/SNPs.
mask = genes["seqname"].astype(str).str.fullmatch(r"chr([1-9]|1[0-9]|2[0-2])")
genes_1_22 = genes[mask].copy()
genes_1_22["seqname"] = genes_1_22["seqname"].str.replace("chr", "", regex=False)

tss_df = genes_1_22[["seqname", "Start", "End", "gene_name"]].copy()
tss_df.columns = ["Chromosome", "Start", "End", "gene"]
tss_gr = pr.PyRanges(tss_df)

# -----------------------------
# 7) Intersect opposite anchors with promoter windows
# -----------------------------
opp_tss_gr = opp_gr.join(tss_gr)
opp_tss = opp_tss_gr.as_df()

# Keep only needed columns (loop_id → gene)
opp_tss = opp_tss[["loop_id", "gene"]].drop_duplicates()

# -----------------------------
# 8) Merge: SNP-hit anchors + (loop_id → gene)
# -----------------------------
hit = snp_hits.merge(opp_tss, on="loop_id", how="inner")

# Save final mapping: GWAS rsID → loop → gene (plus anchor info if desired)
hit.to_csv("risk_loop.csv", index=False)

print(f"Saved {hit.shape[0]} SNP-loop-gene links to risk_loop.csv")






