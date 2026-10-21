
library(data.table)
library(dplyr)
library(coloc)

# --------------------------------------------------
# 1. Set analysis parameters and load GWAS data
# --------------------------------------------------
celltypes <- c("CD4T", "CD8T")

n_case <- 16159
n_control <- 25386
n_eqtl <- 1412

gwas <- data.frame(fread("NG_2021_JT.ma"))

all_result <- list()


# --------------------------------------------------
# 2. Perform colocalization for each T-cell subset
# --------------------------------------------------
for (ct in celltypes) {

    message("Processing: ", ct)

    # Load T1D risk eGenes and eQTL meta-analysis results.
    risk <- read.csv(
        paste0(ct, "_eqtl_meta_results_GWASrisk.csv")
    )

    re_g <- unique(risk$Gene)

    eqtl_0 <- fread(
        paste0(ct, "_eqtl_meta_results.csv")
    )

    eqtl_0$beta <- eqtl_0$Meta_Beta
    eqtl_0$se <- eqtl_0$Meta_SE
    eqtl_0$gene <- eqtl_0$Gene
    eqtl_0$snps <- eqtl_0$rsID

    # Load gene-specific expression standard deviations.
    sdY_df <- read.csv(
        paste0(ct, "_sdY.csv")
    )

    rownames(sdY_df) <- sdY_df$gene

    result_list <- list()
    snp_result_list <- list()


    # --------------------------------------------------
    # 3. Run gene-wise GWAS-eQTL colocalization
    # --------------------------------------------------
    for (gene_i in re_g) {

        eqtl <- eqtl_0[gene == gene_i]

        if (nrow(eqtl) == 0) next

        res <- tryCatch(

            coloc.abf(

                dataset1 = list(
                    beta = gwas$b,
                    varbeta = gwas$se^2,
                    snp = gwas$ID,
                    type = "cc",
                    s = n_case / (n_case + n_control),
                    N = n_case + n_control
                ),

                dataset2 = list(
                    beta = eqtl$beta,
                    varbeta = eqtl$se^2,
                    snp = eqtl$snps,
                    sdY = sdY_df[gene_i, "sdY"],
                    type = "quant",
                    N = n_eqtl
                )

            ),

            error = function(e) {

                message(
                    "Skip ", ct, " | ", gene_i,
                    " | ", e$message
                )

                return(NULL)
            }
        )

        if (is.null(res)) next


        # Extract SNP-level posterior probabilities.
        snp_res <- as.data.frame(res$results)

        top_idx <- which.max(snp_res$SNP.PP.H4)

        top_shared_SNP <- snp_res$snp[top_idx]
        top_SNP_PP_H4 <- snp_res$SNP.PP.H4[top_idx]


        # Identify lead GWAS and eQTL variants.
        lead_GWAS_SNP <- gwas$SNP[
            which.max(abs(gwas$b / gwas$se))
        ]

        lead_eQTL_SNP <- eqtl$snps[
            which.max(abs(eqtl$beta / eqtl$se))
        ]


        # --------------------------------------------------
        # 4. Summarize gene-level colocalization results
        # --------------------------------------------------
        PP3 <- as.numeric(res$summary["PP.H3.abf"])
        PP4 <- as.numeric(res$summary["PP.H4.abf"])

        result_list[[gene_i]] <- data.frame(

            ct = ct,
            gene = gene_i,

            nsnp = as.numeric(
                res$summary["nsnps"]
            ),

            PP0 = as.numeric(
                res$summary["PP.H0.abf"]
            ),

            PP1 = as.numeric(
                res$summary["PP.H1.abf"]
            ),

            PP2 = as.numeric(
                res$summary["PP.H2.abf"]
            ),

            PP3 = PP3,
            PP4 = PP4,

            PP4_conditional = PP4 / (PP3 + PP4),

            lead_GWAS_SNP = lead_GWAS_SNP,
            lead_eQTL_SNP = lead_eQTL_SNP,

            top_shared_SNP = top_shared_SNP,
            top_SNP_PP_H4 = top_SNP_PP_H4,

            status = "success"
        )


        # Store SNP-level colocalization results.
        snp_res$ct <- ct
        snp_res$gene <- gene_i

        snp_result_list[[gene_i]] <- snp_res

    }


    # --------------------------------------------------
    # 5. Export cell-type-specific results
    # --------------------------------------------------

    # Gene-level colocalization summary.
    coloc_result <- bind_rows(result_list) %>%
        arrange(desc(PP4))

    fwrite(
        coloc_result,
        paste0(ct, "_coloc_summary.csv")
    )

    # SNP-level posterior probabilities.
    coloc_snp_result <- bind_rows(snp_result_list)

    fwrite(
        coloc_snp_result,
        paste0(ct, "_coloc_SNP.csv.gz")
    )

    all_result[[ct]] <- coloc_result

}


# --------------------------------------------------
# 6. Combine colocalization results across T-cell subsets
# --------------------------------------------------
all_coloc_result <- bind_rows(all_result) %>%
    arrange(desc(PP4))

fwrite(
    all_coloc_result,
    "all_celltypes_coloc_summary.csv"
)
