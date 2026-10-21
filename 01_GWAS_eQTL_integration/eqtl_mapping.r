library(MatrixEQTL)

# =========================
# 1) Parse command-line arguments
#    arg1: cell type
#    arg2: chromosome
# =========================
args <- commandArgs(trailingOnly = TRUE)
ct  <- args[1]
chr <- args[2]

# Set eQTL model (linear regression)
useModel <- modelLINEAR

# =========================
# 2) Define input/output files
# =========================
# Genotype matrix and SNP positions
SNP_file_name          <- paste0("./demo/onek1k_imputed_merged_QC_filter_012_chr", chr, ".txt")
snps_location_file_name <- "./demo/snploc.txt"

# Gene expression matrix and gene positions
expression_file_name    <- paste0("./demo/GE.", ct, ".txt")
gene_location_file_name <- "./demo/geneloc.txt"

# Covariates matrix (set to character() if no covariates are used)
covariates_file_name <- "./demo/Covariates_cli.txt"

# Temporary output files (MatrixEQTL writes results to disk and reads them back)
output_file_name_cis <- tempfile()
output_file_name_tra <- tempfile()

# Save only associations below these p-value thresholds
pvOutputThreshold_cis <- 2e-2
pvOutputThreshold_tra <- 1e-2

# Error covariance matrix (numeric() corresponds to identity)
errorCovariance <- numeric()

# cis distance threshold: ±1 Mb (1,000,000 bp)
cisDist <- 1e6

# =========================
# 3) Load genotype data (SNPs × samples)
# =========================
snps <- SlicedData$new()
snps$fileDelimiter      <- " "   # Use "\t" for tab-delimited files; use "" for any whitespace
snps$fileOmitCharacters <- "NA"  # Missing value token
snps$fileSkipRows       <- 1     # Skip header row
snps$fileSkipColumns    <- 1     # Skip row-name column (SNP IDs)
snps$fileSliceSize      <- 2000  # Read data in chunks
snps$LoadFile(SNP_file_name)

# =========================
# 4) Load gene expression data (genes × samples)
# =========================
gene <- SlicedData$new()
gene$fileDelimiter      <- " "
gene$fileOmitCharacters <- "NA"
gene$fileSkipRows       <- 1
gene$fileSkipColumns    <- 1     # Skip row-name column (gene IDs)
gene$fileSliceSize      <- 2000
gene$LoadFile(expression_file_name)

# =========================
# 5) Load covariates (covariates × samples)
# =========================
cvrt <- SlicedData$new()
cvrt$fileDelimiter      <- " "
cvrt$fileOmitCharacters <- "NA"
cvrt$fileSkipRows       <- 1
cvrt$fileSkipColumns    <- 1     # Skip row-name column (covariate IDs)
if (length(covariates_file_name) > 0) {
  cvrt$LoadFile(covariates_file_name)
}

# =========================
# 6) Load SNP and gene genomic positions
# =========================
snpspos <- read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
snpspos$snpid <- rownames(snpspos)
snpspos <- snpspos[, c("snpid", "chr", "pos")]

genepos <- read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos$geneid <- rownames(genepos)
genepos <- genepos[, c("geneid", "chr", "left", "right")]

# =========================
# 7) Match sample IDs across genotype, expression, and covariates
# =========================
snp_samples  <- snps$columnNames
gene_samples <- gene$columnNames
common_samples <- intersect(snp_samples, gene_samples)

# Subset genotype and expression matrices to shared samples (same order)
snps$ColumnSubsample(match(common_samples, snps$columnNames))
gene$ColumnSubsample(match(common_samples, gene$columnNames))

# Subset covariates to shared samples (if covariates are provided)
cvrt$ColumnSubsample(match(common_samples, cvrt$columnNames))

# =========================
# 8) Run MatrixEQTL
# =========================
me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE
)

# Remove temporary result files
unlink(output_file_name_tra)
unlink(output_file_name_cis)

# =========================
# 9) Save cis-eQTL results
# =========================
res <- me$cis$eqtls
write.csv(res, paste0(ct, "_chr", chr, "_cis_eQTL.csv"), row.names = FALSE)