#!/usr/bin/env bash

# Activate conda environment (adjust if your setup uses conda.sh)
source activate eqtl_env

for chr in {1..22} X; do
  echo "Processing chromosome chr${chr}..."
  Rscript run_eQTL.r CD4T "${chr}"
done