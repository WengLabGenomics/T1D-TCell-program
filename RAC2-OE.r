library(data.table)
library(stringr)
library(edgeR)

# -----------------------------
# 1) Load raw count matrix (featureCounts-style)
# -----------------------------
df <- as.data.frame(fread("overexpression_counts.txt"))

# Use Geneid as rownames and drop non-count annotation columns
rownames(df) <- df$Geneid
df <- df[, 7:ncol(df)]  # keep only count columns (adjust if your file has a different layout)

# Rename samples: 4 WT and 4 RAC2
colnames(df) <- c(paste0("WT-", 1:4), paste0("RAC2-", 1:4))

# -----------------------------
# 2) Build sample metadata
# -----------------------------
meta <- data.frame(row.names = colnames(df))
meta$group <- str_split_fixed(colnames(df), "-", 2)[, 1]

# Ensure WT is the reference group
meta$group <- factor(meta$group, levels = c("WT", "RAC2"))
meta$group <- relevel(meta$group, ref = "WT")

# -----------------------------
# 3) Create DGEList and filter lowly expressed genes
# -----------------------------
dge <- DGEList(counts = df, group = meta$group)

# Keep genes with CPM > 1 in at least 2 samples
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Recompute library sizes after filtering
dge$samples$lib.size <- colSums(dge$counts)

# TMM normalization
dge <- calcNormFactors(dge)

# -----------------------------
# 4) Design matrix (no intercept): WT and RAC2 columns
# -----------------------------
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)

# -----------------------------
# 5) Estimate dispersions and fit NB GLM
# -----------------------------
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)

# Test RAC2 vs WT (RAC2 - WT)
lrt <- glmLRT(fit, contrast = c(-1, 1))

# -----------------------------
# 6) Export differential expression results
# -----------------------------
res <- topTags(lrt, n = nrow(dge))$table
write.csv(res, "RAC2_DEG.csv", row.names = TRUE)






