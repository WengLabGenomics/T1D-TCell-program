# ============================================================
# eGene program scoring across cohorts
# 1) Cross-sectional scRNA-seq (pseudo-bulk from Seurat object)
# 2) Longitudinal microarray (RMA-normalized expression)
# 3) Teplizumab cohort (bulk RNA-seq counts)
# ============================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(stringr)
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(Matrix)
})

# ----------------------------
# Helper: compute eGene program score from an expression matrix
# - expr_mat: genes x samples (numeric)
# - genes: character vector of eGenes
# ----------------------------
compute_program_score <- function(expr_mat, genes, logCPM = TRUE) {
  genes_use <- intersect(genes, rownames(expr_mat))
  if (length(genes_use) == 0) stop("No overlap between eGenes and expression matrix rownames.")

  if (logCPM) {
    expr_mat <- edgeR::cpm(expr_mat, log = TRUE, prior.count = 1)
  }
  X <- t(expr_mat[genes_use, , drop = FALSE])  # samples x genes
  rowMeans(X, na.rm = TRUE)
}

# ============================================================
# 0) Load eGene list
# ============================================================
eg_df <- read.csv("eGene-list.csv", row.names = 1)   # adjust filename if needed
eGenes <- eg_df[[1]]                                # assumes first column is the gene symbol list

# ============================================================
# 1) Cross-sectional scRNA-seq cohort (46 T1D vs 31 HC)
#    Build pseudo-bulk counts by summing UMI counts per sample
# ============================================================
obj <- readRDS("T1D_Seurat_Object_Final.rds")

# Keep only T cells (CD4_T and CD8_T)
obj_t <- subset(obj, subset = Cluster_Annotation_Merged %in% c("CD4_T", "CD8_T"))

# Extract raw counts (genes x cells)
counts <- GetAssayData(obj_t, assay = "RNA", slot = "counts")

# Define the grouping variable used to aggregate cells into samples
group_col <- "Sample_ID"
samples <- obj_t[[group_col]][, 1]

# Pseudo-bulk: sum counts across cells for each sample
pb_counts <- sapply(split(seq_len(ncol(counts)), samples), function(idx) {
  Matrix::rowSums(counts[, idx, drop = FALSE])
})
pb_counts <- as.matrix(pb_counts)  # genes x samples

# (Optional) save pseudo-bulk counts
write.csv(pb_counts, "./T_cell_pseudobulk_counts.csv")

# Compute program score (logCPM from pseudo-bulk counts)
score_sc <- compute_program_score(pb_counts, eGenes, logCPM = TRUE)

# Build metadata (HC vs T1D). Adjust parsing to your sample naming scheme.
meta_sc <- data.frame(sample = colnames(pb_counts), score = score_sc, stringsAsFactors = FALSE)
# Example: if sample IDs are like "H.xxx" or "T1D.xxx"
meta_sc$group <- str_split_fixed(meta_sc$sample, "\\.", 2)[, 1]
meta_sc$group <- recode(meta_sc$group, "H" = "HC")
meta_sc$group <- factor(meta_sc$group, levels = c("HC", "T1D"))


# Plot
p_sc <- ggplot(meta_sc, aes(x = group, y = score, color = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, color = "black") +
  theme_classic() +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 0)) +
  labs(x = NULL, y = "eGene program score") +
  stat_compare_means(comparisons = list(c("HC", "T1D")),
                     method = "t.test", label = "p.signif")
print(p_sc)

# ============================================================
# 2) Longitudinal microarray cohort (GSE43488)
#    RMA normalize CEL files, map probes -> gene symbols,
#    compute program score, and compare groups by t-tests
# ============================================================
suppressPackageStartupMessages({
  library(oligo)
  library(data.table)
})

celfiles <- list.files("GSE43488", full.names = TRUE)
rawData <- read.celfiles(celfiles)
eset <- rma(rawData)
expr_matrix <- exprs(eset)  # probes x samples

# Load platform annotation table (must contain probe ID and gene symbol columns)
an <- as.data.frame(fread("./GPL13667-15572.txt", header = TRUE))
an <- an[an$Gene.Symbol != "---", ]
rownames(an) <- an$ID

keep <- intersect(rownames(expr_matrix), an$ID)
expr_matrix2 <- expr_matrix[keep, , drop = FALSE]
rownames(expr_matrix2) <- an[keep, "Gene.Symbol"]

# Collapse duplicated gene symbols by keeping the first occurrence
expr_matrix3 <- expr_matrix2[!duplicated(rownames(expr_matrix2)), , drop = FALSE]

# (Optional) save normalized expression matrix
write.csv(expr_matrix3, "./high-risk_normalized_microarray.csv")

# Sample annotation
sa <- data.frame(sample = colnames(expr_matrix3), stringsAsFactors = FALSE)

# Example parsing: adjust to your sample naming
tmp <- str_split_fixed(sa$sample, "\\.", 3)[, 1]
case <- str_split_fixed(tmp, "_", 3)[, 2]
pro  <- str_split_fixed(tmp, "_", 3)[, 3]
sa$sample_ID <- paste0(case, "_", pro)
rownames(sa) <- sa$sample_ID

metadata <- read.csv("DIPP_meta.csv", row.names = 1)
sa$Time_from_dx_months <- metadata[rownames(sa), "Time.from.T1D.diagnosis..months."]

sa$group <- NA_character_
sa[grepl("HC", sa$sample_ID), "group"] <- "Healthy control"
sa[!is.na(sa$Time_from_dx_months) &
     sa$Time_from_dx_months > -24 & sa$Time_from_dx_months < 0, "group"] <- 'The two years before the diagnosis of T1D'
sa[!is.na(sa$Time_from_dx_months) &
     sa$Time_from_dx_months >= 0, "group"] <- 'At the time of clinical diagnosis of T1D'

# Compute program score (microarray already normalized; do NOT logCPM)
score_ma <- compute_program_score(expr_matrix3, eGenes, logCPM = FALSE)

sd_ma <- data.frame(
  sample_ID = colnames(expr_matrix3),
  score = score_ma,
  stringsAsFactors = FALSE
)
# Align sample IDs
sd_ma$sample_ID <- str_split_fixed(sd_ma$sample_ID, "\\.", 3)[, 1]
sd_ma$sample_ID <- paste0(str_split_fixed(sd_ma$sample_ID, "_", 3)[, 2], "_",
                          str_split_fixed(sd_ma$sample_ID, "_", 3)[, 3])
rownames(sd_ma) <- sd_ma$sample_ID
sd_ma <- cbind(sd_ma, sa[rownames(sd_ma), , drop = FALSE])
sd_ma <- sd_ma[!is.na(sd_ma$group), , drop = FALSE]

sd_ma$group <- factor(
  sd_ma$group,
  levels = c("Healthy control", 'The two years before the diagnosis of T1D', 'At the time of clinical diagnosis of T1D')
)

# Plot
p_ma <- ggplot(sd_ma, aes(x = group, y = score, color = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, color = "black") +
  theme_classic() +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(x = NULL, y = "eGene program score") +
  stat_compare_means(
    comparisons = list(
      c("Healthy control", 'The two years before the diagnosis of T1D'),
      c( 'The two years before the diagnosis of T1D','At the time of clinical diagnosis of T1D'),
      c("Healthy control", 'At the time of clinical diagnosis of T1D')
    ),
    method = "t.test",
    label = "p.signif"
  )
print(p_ma)

# ============================================================
# 3) Teplizumab cohort (GSE85531)
#    Compute program score from raw counts using logCPM,
#    then test longitudinal change with a mixed model
# ============================================================
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

df_counts <- read.table(
  "GSE85531_combined_ABATE_raw_counts_after_qc_063016_not_collapsed.txt",
  header = TRUE, row.names = 1, check.names = FALSE
)

# Map Ensembl IDs -> gene symbols
gene_anno <- as.data.frame(fread("Human.GRCh38.p13.annot.tsv"))
gene_anno <- gene_anno[gene_anno$EnsemblGeneID %in% rownames(df_counts), ]
gene_anno <- gene_anno[!duplicated(gene_anno$EnsemblGeneID), ]
rownames(gene_anno) <- gene_anno$EnsemblGeneID

df_counts <- df_counts[rownames(gene_anno), , drop = FALSE]
rownames(df_counts) <- gene_anno$Symbol

meta_tp <- read.csv("SraRunTable.csv", stringsAsFactors = FALSE)
rownames(meta_tp) <- meta_tp$samplelabel

# Keep responders/non-responders only
meta_tp2 <- meta_tp[meta_tp$status %in% c("R", "NR"), , drop = FALSE]
df_tp <- df_counts[, rownames(meta_tp2), drop = FALSE]

# logCPM transform
logcpm_tp <- edgeR::cpm(df_tp, log = TRUE, prior.count = 1)
score_tp <- compute_program_score(logcpm_tp, eGenes, logCPM = FALSE)

sd_tp <- data.frame(
  samplelabel = colnames(df_tp),
  score = score_tp,
  stringsAsFactors = FALSE
)
sd_tp <- cbind(sd_tp, meta_tp2[sd_tp$samplelabel, , drop = FALSE])

# Keep paired subjects with both baseline (0) and month 6
sd_tp <- sd_tp[sd_tp$visitmonth %in% c(0, 6), , drop = FALSE]
paired_ids <- intersect(sd_tp$itn_name[sd_tp$visitmonth == 0],
                        sd_tp$itn_name[sd_tp$visitmonth == 6])
sd_tp <- sd_tp[sd_tp$itn_name %in% paired_ids, , drop = FALSE]

sd_tp$time <- factor(sd_tp$visitmonth, levels = c(0, 6), labels = c("0 month", "6 month"))
sd_tp$status <- factor(sd_tp$status, levels = c("NR", "R"))

# Plot with paired lines and boxplots
p_tp <- ggplot(sd_tp, aes(x = time, y = score, color = status)) +
  geom_line(aes(group = itn_name), alpha = 0.5, linewidth = 0.6, color = "gray55") +
  geom_point(size = 1.8, color = "black") +
  geom_boxplot(aes(group = time), fill = "white", width = 0.4, outlier.shape = NA, alpha = 0.6) +
  facet_wrap(~ status) +
  theme_classic(base_size = 13) +
  labs(x = "Visit month", y = "eGene program score (logCPM)") +
  theme(panel.grid = element_blank())
print(p_tp)

# Mixed model:
# score ~ time * status + (1 | subject)
fit_noint <- lmer(score ~ time + status + (1 | itn_name), data = sd_tp, REML = FALSE)
fit_full  <- lmer(score ~ time * status + (1 | itn_name), data = sd_tp, REML = FALSE)

# Likelihood-ratio test for interaction (time-by-response)
print(anova(fit_noint, fit_full))






