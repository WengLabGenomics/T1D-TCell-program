library(data.table)
library(dplyr)
library(ggplot2)

# --------------------------------------------------
# 1. Load Hi-C, eQTL, and ATAC-seq results
# --------------------------------------------------
ct <- "CD4T"  # or "CD8T"

cell_label <- if (ct == "CD4T") "CD4+ T cells" else "CD8+ T cells"

loop <- fread(
    paste0(ct, "_mustache_SNP_to_gene_candidates_2kb_gene.csv")
)

meta <- fread(
    paste0(ct, "_eqtl_meta_results_GWASrisk.csv")
)

ATAC <- fread(
    paste0("SNP_to_ATAC_bw_overlap_", ct, ".csv")
)


# --------------------------------------------------
# 2. Identify risk-associated Hi-C loops
# --------------------------------------------------

# Retain loops linked to T1D risk eGenes.
loop2 <- loop %>%
    filter(gene %in% meta$Gene)

# Retain loops containing SNPs with ATAC-seq signals
# at or above the 75th percentile (excluding zero signal).
loop3 <- loop2 %>%
    filter(rsID %in% ATAC$rsID[ATAC$open_P75 == 1])

# Export filtered SNP-gene links supported by Hi-C and ATAC-seq.
write.csv(
    loop3,
    paste0(ct, "_risk_loop.csv"),
    row.names = FALSE
)


# --------------------------------------------------
# 3. Prepare genomic coordinates for plotting
# --------------------------------------------------

# Retain one representative SNP-overlapping anchor per loop.
loops <- loop3 %>%
    distinct(loop_id, .keep_all = TRUE) %>%
    mutate(
        chrom = sub("^chr", "", as.character(chr)),
        position = (as.numeric(start) + as.numeric(end)) / 2,
        fdr = as.numeric(FDR)
    ) %>%
    filter(
        chrom %in% as.character(1:22),
        !is.na(fdr),
        fdr >= 0
    )


# --------------------------------------------------
# 4. Define chromosome lengths (hg38)
# --------------------------------------------------
chrom_sizes <- data.frame(
    chrom = as.character(1:22),
    size = c(
        248956422, 242193529, 198295559, 190214555,
        181538259, 170805979, 159345973, 145138636,
        138394717, 133797422, 135086622, 133275309,
        114364328, 107043718, 101991189, 90338345,
        83257441, 80373285, 58617616, 64444167,
        46709983, 50818468
    )
)

# Calculate cumulative genomic positions across chromosomes.
chrom_sizes <- chrom_sizes %>%
    mutate(
        chr_start = cumsum(size) - size,
        chr_end = chr_start + size,
        chr_center = chr_start + size / 2
    )


# --------------------------------------------------
# 5. Calculate genome-wide plotting coordinates
# --------------------------------------------------
loops <- loops %>%
    left_join(chrom_sizes, by = "chrom") %>%
    mutate(
        pos_cum = position + chr_start,
        neg_log10_fdr = -log10(
            pmax(fdr, .Machine$double.xmin)
        )
    )


# --------------------------------------------------
# 6. Plot genomic distribution of T1D-associated Hi-C loops
# --------------------------------------------------
p <- ggplot(
    loops,
    aes(x = pos_cum, y = neg_log10_fdr)
) +

    geom_point(
        color = "#8694CA",
        size = 2,
        alpha = 1
    ) +

    # Display all 22 autosomes.
    scale_x_continuous(
        breaks = chrom_sizes$chr_center,
        labels = chrom_sizes$chrom,
        limits = c(0, max(chrom_sizes$chr_end)),
        expand = c(0, 0)
    ) +

    # FDR significance threshold.
    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed"
    ) +

    labs(
        x = "Chromosome",
        y = expression(-log[10](FDR)),
        title = paste0(
            "Genomic distribution of T1D-associated Hi-C loops in ",
            cell_label
        )
    ) +

    theme_classic(base_size = 13) +

    theme(
        aspect.ratio = 0.2,
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = 0,
            hjust = 0.5
        ),
        plot.title = element_text(
            hjust = 0.5
        )
    )

print(p)

