# =========================
# Fixed-effect meta-analysis for cis-eQTL summary statistics
# - Merge per-study eQTL tables by SNP–gene pair
# - Keep eQTLs supported by at least 2 studies (i.e., <3 missing betas out of 4)
# - Run fixed-effect meta-analysis using metafor::rma()
# - Adjust meta P-values using BH FDR
# =========================

library(data.table)
library(dplyr)
library(metafor)

# -------------------------
# 1) Load per-study eQTL summary statistics
#    Required columns (per file): rsID, Gene, beta, se, p
# -------------------------

# pathogen_NC_2022
df1 <- as.data.frame(fread("pathogen_NC_2022_CD4T.csv"))
df1$eqtl <- paste0(df1$rsID, "_", df1$Gene)  # Define unique SNP–gene key

# Nature_2023
df2 <- as.data.frame(fread("Nature_2023_CD4T.csv"))
df2$eqtl <- paste0(df2$rsID, "_", df2$Gene)

# onek1k
df3 <- as.data.frame(fread("onek1k_CD4T.csv"))
df3$eqtl <- paste0(df3$rsID, "_", df3$Gene)

# Viral_Science_2021
df4 <- as.data.frame(fread("Viral_Science_2021_CD4T.csv"))
df4$eqtl <- paste0(df4$rsID, "_", df4$Gene)

# -------------------------
# 2) Harmonize and rename study-specific effect columns
#    Keep SNP identifiers from the first study (if available)
# -------------------------
paper1 <- df1 %>%
  dplyr::select(eqtl, rsID, SNP1, SNP2, Gene, beta, se, p) %>%
  rename_with(~ paste0(., "_p1"), .cols = c(beta, se, p))

paper2 <- df2 %>%
  dplyr::select(eqtl, beta, se, p) %>%
  rename_with(~ paste0(., "_p2"), .cols = c(beta, se, p))

paper3 <- df3 %>%
  dplyr::select(eqtl, beta, se, p) %>%
  rename_with(~ paste0(., "_p3"), .cols = c(beta, se, p))

paper4 <- df4 %>%
  dplyr::select(eqtl, beta, se, p) %>%
  rename_with(~ paste0(., "_p4"), .cols = c(beta, se, p))

# -------------------------
# 3) Merge all studies by the SNP–gene key
# -------------------------
merged <- paper1 %>%
  full_join(paper2, by = "eqtl") %>%
  full_join(paper3, by = "eqtl") %>%
  full_join(paper4, by = "eqtl")

# Remove duplicated keys (keep the first occurrence)
merged_clean <- merged[!duplicated(merged$eqtl), ]

# Keep eQTLs observed in at least 2 studies
# (i.e., fewer than 3 missing beta values out of 4)
merged_clean <- merged_clean[
  apply(
    merged_clean[, c("beta_p1", "beta_p2", "beta_p3", "beta_p4")],
    1,
    function(x) sum(is.na(x)) < 3
  ),
]

# -------------------------
# 4) Meta-analysis helper (fixed-effect model)
#    - Uses metafor::rma(method = "FE")
#    - Skips missing beta/se pairs automatically
# -------------------------
run_meta <- function(betas, ses) {
  # Identify indices where both beta and se are available
  valid_idx <- which(!is.na(betas) & !is.na(ses))

  # Optionally enforce a minimum number of studies
  # (Uncomment if you want >=3 studies)
  # if (length(valid_idx) < 3) return(c(NA_real_, NA_real_, NA_real_))

  betas <- betas[valid_idx]
  ses   <- ses[valid_idx]

  # Run fixed-effect meta-analysis; return NA if it fails
  res <- tryCatch(
    rma(yi = betas, sei = ses, method = "FE"),
    error = function(e) NULL
  )

  if (is.null(res)) {
    return(c(NA_real_, NA_real_, NA_real_))
  } else {
    return(c(res$beta, res$se, res$pval))
  }
}

# -------------------------
# 5) Run meta-analysis for each eQTL (row-wise)
# -------------------------
meta_results <- t(apply(merged_clean, 1, function(row) {
  betas <- as.numeric(row[c("beta_p1", "beta_p2", "beta_p3", "beta_p4")])
  ses   <- as.numeric(row[c("se_p1",   "se_p2",   "se_p3",   "se_p4")])
  run_meta(betas, ses)
}))

meta_results_df <- as.data.frame(meta_results)
colnames(meta_results_df) <- c("Meta_Beta", "Meta_SE", "Meta_Pval")

# -------------------------
# 6) Combine results and adjust for multiple testing
# -------------------------
final <- cbind(
  merged_clean,
  meta_results_df
)

final$fdr <- p.adjust(final$Meta_Pval, method = "fdr")

# -------------------------
# 7) Write output
# -------------------------
write.csv(final, "CD4T_eqtl_meta_results.csv", row.names = FALSE)