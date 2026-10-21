import numpy as np
import pandas as pd
import pyranges as pr
import pyBigWig


# --------------------------------------------------
# 1. Configure cell type and analysis thresholds
# --------------------------------------------------
ct = "CD4"  # Change to "CD8" to analyze CD8+ T cells.

GWAS_P_THRESHOLD = 5e-8
MIN_DIST = 20_000
MAX_DIST = 5_000_000
PROMOTER_FLANK = 2_000


def normalize_chromosome(series):
    """Convert chromosome labels such as 'chr1' and 1 to '1'."""
    return (
        series.astype(str)
        .str.replace(r"^chr", "", regex=True)
        .str.replace(r"\.0$", "", regex=True)
    )


# --------------------------------------------------
# 2. Load Hi-C loops and construct loop identifiers
# --------------------------------------------------
loops_raw = pd.read_table(f"mustache_SCALE_{ct}.tsv")

# Accommodate alternative chromosome column names in Mustache output.
loops_raw = loops_raw.rename(
    columns={
        "BIN1_CHROMOSOME": "BIN1_CHR",
        "BIN2_CHR": "BIN2_CHROMOSOME",
    }
)

loop_columns = {
    "BIN1_CHR": "chrA",
    "BIN1_START": "A_start",
    "BIN1_END": "A_end",
    "BIN2_CHROMOSOME": "chrB",
    "BIN2_START": "B_start",
    "BIN2_END": "B_end",
    "FDR": "FDR",
    "DETECTION_SCALE": "DETECTION_SCALE",
}

loops = loops_raw[list(loop_columns)].rename(columns=loop_columns).copy()

loops["chrA"] = normalize_chromosome(loops["chrA"])
loops["chrB"] = normalize_chromosome(loops["chrB"])

# Use both genomic anchors to define a unique loop identifier.
loop_coordinates = [
    "chrA", "A_start", "A_end",
    "chrB", "B_start", "B_end"
]

loops["loop_id"] = (
    loops[loop_coordinates].astype(str).agg("_".join, axis=1)
)

# Calculate the distance between loop anchor midpoints.
midA = (loops["A_start"] + loops["A_end"]) // 2
midB = (loops["B_start"] + loops["B_end"]) // 2

same_chr = loops["chrA"] == loops["chrB"]

loops["distance"] = (
    (midA - midB).abs().where(same_chr, np.nan)
)

# Apply distance limits to cis loops.
# Retain trans loops as in the original analysis.
loops_f = loops.loc[
    (~same_chr) |
    loops["distance"].between(MIN_DIST, MAX_DIST)
].copy()


# --------------------------------------------------
# 3. Convert both loop anchors to genomic intervals
# --------------------------------------------------

# Extract the first loop anchor.
A = loops_f[["chrA", "A_start", "A_end", "loop_id"]].copy()
A.columns = ["chr", "start", "end", "loop_id"]
A["side"] = "A"

# Extract the second loop anchor.
B = loops_f[["chrB", "B_start", "B_end", "loop_id"]].copy()
B.columns = ["chr", "start", "end", "loop_id"]
B["side"] = "B"

# Combine both anchors for genomic overlap analysis.
anchors = pd.concat([A, B], ignore_index=True)

anchors_df = anchors.rename(
    columns={
        "chr": "Chromosome",
        "start": "Start",
        "end": "End"
    }
)

anchors_gr = pr.PyRanges(anchors_df)


# --------------------------------------------------
# 4. Load GWAS-significant variants with positive beta
# --------------------------------------------------
ma = pd.read_csv("NG_2021_JT.ma")

# Select genome-wide significant variants with positive effect sizes.
# Interpret beta > 0 as risk-increasing only after checking allele orientation.
ma = ma.loc[
    (ma["p_value"] < GWAS_P_THRESHOLD) &
    (ma["beta"] > 0)
].copy()

# Extract genomic coordinates and SNP identifiers.
snps_df = ma[
    ["chromosome", "base_pair_location", "V3"]
].rename(
    columns={
        "chromosome": "Chromosome",
        "base_pair_location": "Position",
        "V3": "rsID"
    }
).copy()

# Standardize chromosome names and retain autosomes.
snps_df["Chromosome"] = normalize_chromosome(
    snps_df["Chromosome"]
)

snps_df = snps_df.loc[
    snps_df["Chromosome"].isin(
        [str(i) for i in range(1, 23)]
    )
].copy()

# Assume GWAS positions are 1-based.
# Convert SNP positions to 0-based, half-open intervals.
snps_df["Start"] = snps_df["Position"].astype(int) - 1
snps_df["End"] = snps_df["Start"] + 1

snps_df = snps_df[
    ["Chromosome", "Start", "End", "rsID"]
].copy()

snps_gr = pr.PyRanges(snps_df)


# --------------------------------------------------
# 5. Identify GWAS variants overlapping Hi-C loop anchors
# --------------------------------------------------

# Find SNPs overlapping either anchor of a Hi-C loop.
snp_hits = anchors_gr.join(snps_gr).as_df()

# Standardize column names for downstream analysis.
snp_hits = snp_hits.rename(
    columns={
        "Chromosome": "chr",
        "Start": "start",
        "End": "end",
        "Start_b": "snp_start",
        "End_b": "snp_end"
    }
)[
    [
        "chr", "start", "end",
        "loop_id", "side",
        "snp_start", "snp_end", "rsID"
    ]
]


# --------------------------------------------------
# 6. Identify the opposite anchor for each SNP hit
# --------------------------------------------------

# Construct a lookup table for the opposite loop anchor.
opposite_lookup = anchors.rename(
    columns={
        "chr": "opp_chr",
        "start": "opp_start",
        "end": "opp_end",
        "side": "opp_side"
    }
).copy()

# A SNP at anchor A is linked to anchor B, and vice versa.
opposite_lookup["side"] = (
    opposite_lookup["opp_side"].map(
        {"A": "B", "B": "A"}
    )
)

opposite_lookup = opposite_lookup.drop_duplicates(
    ["loop_id", "side"]
)

# Assign the opposite anchor to each SNP-overlapping anchor.
snp_hits = snp_hits.merge(
    opposite_lookup[
        [
            "loop_id", "side",
            "opp_chr", "opp_start",
            "opp_end", "opp_side"
        ]
    ],
    on=["loop_id", "side"],
    how="left",
    validate="many_to_one"
)

# Prepare opposite-anchor intervals for promoter overlap analysis.
opp_df = snp_hits[
    [
        "opp_chr", "opp_start", "opp_end",
        "loop_id", "opp_side"
    ]
].dropna().drop_duplicates().rename(
    columns={
        "opp_chr": "Chromosome",
        "opp_start": "Start",
        "opp_end": "End"
    }
).copy()

opp_df[["Start", "End"]] = (
    opp_df[["Start", "End"]].astype(int)
)

opp_gr = pr.PyRanges(opp_df)


# --------------------------------------------------
# 7. Map opposite loop anchors to gene promoters (+/- 2 kb)
# --------------------------------------------------

# Load gene annotations containing transcription start sites.
gtf_genes = pd.read_csv("gencode_gtf.csv")

# Standardize chromosome names and retain autosomal genes.
gtf_genes["Chromosome"] = normalize_chromosome(
    gtf_genes["seqname"]
)

gtf_genes = gtf_genes.loc[
    gtf_genes["Chromosome"].isin(
        [str(i) for i in range(1, 23)]
    )
].copy()

# Define promoter regions within 2 kb of each TSS.
# Assume TSS_start is 1-based and convert to 0-based intervals.
tss_pos = gtf_genes["TSS_start"].astype(int)

gtf_genes["Start"] = (
    tss_pos - 1 - PROMOTER_FLANK
).clip(lower=0)

gtf_genes["End"] = (
    tss_pos + PROMOTER_FLANK
)

# Convert promoter coordinates to PyRanges format.
tss_df = gtf_genes[
    ["Chromosome", "Start", "End", "gene_name"]
].rename(
    columns={"gene_name": "gene"}
)

tss_gr = pr.PyRanges(tss_df)

# Identify gene promoters overlapping opposite loop anchors.
opp_tss = opp_gr.join(tss_gr).as_df()

# Match on both loop ID and opposite-anchor side.
# This prevents incorrect SNP-gene assignments within the same loop.
hit = snp_hits.merge(
    opp_tss[
        ["loop_id", "opp_side", "gene"]
    ].drop_duplicates(),
    on=["loop_id", "opp_side"],
    how="inner"
)

# Attach loop significance and genomic distance information.
hit = hit.merge(
    loops_f[
        ["loop_id", "FDR", "distance"]
    ]
    .sort_values("FDR")
    .drop_duplicates("loop_id"),
    on="loop_id",
    how="left"
)

# Rank candidate SNP-gene links by loop FDR and distance.
hit = hit.sort_values(
    ["FDR", "distance"],
    ascending=[True, True]
)

# Export candidate SNP-gene links supported by Hi-C interactions.
hit.to_csv(
    f"{ct}T_mustache_SNP_to_gene_candidates_2kb_gene.csv",
    index=False
)


# --------------------------------------------------
# 8. Quantify ATAC-seq bigWig signal at GWAS variants
# --------------------------------------------------

# Load the merged ATAC-seq bigWig file for the selected T-cell type.
bw = pyBigWig.open(f"merged_{ct}.bw")

if bw is None:
    raise OSError(f"Could not open merged_{ct}.bw")

bw_chromosomes = bw.chroms()


def summarize(chromosome, start, end, mode="mean"):
    """Summarize bigWig signal in a 0-based, half-open genomic interval."""

    chrom = str(chromosome)

    # Match chromosome names to the bigWig reference.
    chrom = next(
        (
            name
            for name in (chrom, f"chr{chrom}")
            if name in bw_chromosomes
        ),
        None
    )

    if chrom is None:
        raise ValueError(
            f"Chromosome {chromosome} is absent from the bigWig file"
        )

    # Retrieve ATAC-seq signal within the SNP interval.
    values = np.asarray(
        bw.values(chrom, int(start), int(end)),
        dtype=float
    )

    # Remove missing or non-finite signal values.
    values = values[np.isfinite(values)]

    if values.size == 0:
        return 0.0

    if mode == "mean":
        return float(values.mean())

    if mode == "max":
        return float(values.max())

    if mode == "sum":
        return float(values.sum())

    raise ValueError(
        "mode must be 'mean', 'max', or 'sum'"
    )


# Prepare GWAS variants for ATAC-seq signal extraction.
atac_snps = snps_df.rename(
    columns={
        "Chromosome": "chromosome",
        "Start": "start",
        "End": "end"
    }
).copy()

# Calculate the mean ATAC-seq signal at each GWAS variant.
try:
    atac_snps["signal"] = [
        summarize(chrom, start, end, mode="mean")
        for chrom, start, end in atac_snps[
            ["chromosome", "start", "end"]
        ].itertuples(index=False, name=None)
    ]

finally:
    bw.close()


# --------------------------------------------------
# 9. Classify variants based on ATAC-seq signal percentiles
# --------------------------------------------------

# Calculate the 65th and 75th percentiles of SNP-level ATAC signal.
thr_65 = atac_snps["signal"].quantile(0.65)
thr_75 = atac_snps["signal"].quantile(0.75)

# Flag variants with nonzero signal above each percentile threshold.
atac_snps["open_P65"] = (
    (atac_snps["signal"] > 0) &
    (atac_snps["signal"] >= thr_65)
).astype(int)

atac_snps["open_P75"] = (
    (atac_snps["signal"] > 0) &
    (atac_snps["signal"] >= thr_75)
).astype(int)

# Export SNP-level ATAC-seq signal and percentile classifications.
atac_snps.to_csv(
    f"SNP_to_ATAC_bw_overlap_{ct}T.csv",
    index=False
)
