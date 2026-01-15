# eQTL module (cis-eQTL mapping + meta-analysis)

This folder contains code to (i) run cis-eQTL mapping per cell type and chromosome using MatrixEQTL and (ii) perform fixed-effect meta-analysis across studies.

## 1) System requirements
- OS: Linux
- R: >= 4.1
- Optional: bash (to run `run_eqtl_mapping.sh`)

## 2) Software dependencies (with versions)
R packages:
- MatrixEQTL (2.3)
- data.table (1.14.8)
- dplyr (1.1.2)
- metafor (4.2_0)

## 3) Versions tested on
Tested with: R 4.1.2 on Linux.

## 4) Non-standard hardware
None required for the demo dataset. (Full-scale analyses may require HPC/large memory depending on input size.)

## 5) Installation + demo (minimal)
Install required R packages: `install.packages(c("data.table","dplyr","metafor")); if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install("MatrixEQTL")`.
A small demo dataset is provided in `demo/` (genotypes: `onek1k_imputed_merged_QC_filter_012_chr*.txt`; SNP positions: `snploc.txt`; expression: `GE.<celltype>.txt`; gene positions: `geneloc.txt`; covariates: `Covariates_cli.txt`).
To run the demo cis-eQTL mapping from `eqtl/`: `Rscript eqtl_mapping.r CD4T 1`, which generates `CD4T_chr1_cis_eQTL.csv` (expected runtime ~1–5 min on a normal desktop).
