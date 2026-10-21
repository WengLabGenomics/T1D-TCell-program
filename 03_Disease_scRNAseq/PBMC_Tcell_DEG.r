library(Seurat)
library(ggplot2)

# --------------------------------------------------
# 1. Load PBMC Seurat object
# --------------------------------------------------
pbmc <- readRDS("T1D_Seurat_Object_Final.rds")

# Use RNA expression data for visualization and DE analysis.
DefaultAssay(pbmc) <- "RNA"

# Join split RNA layers if present (Seurat v5).
if (inherits(pbmc[["RNA"]], "Assay5")) {
    layer_names <- Layers(pbmc[["RNA"]])

    if (any(grepl("^(counts|data)\\.", layer_names))) {
        pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
    }
}


# --------------------------------------------------
# 2. Normalize data and perform dimensionality reduction
# --------------------------------------------------
pbmc <- NormalizeData(pbmc, verbose = FALSE)

pbmc <- FindVariableFeatures(
    pbmc,
    selection.method = "vst",
    nfeatures = 2000
)

pbmc <- ScaleData(
    pbmc,
    features = VariableFeatures(pbmc),
    verbose = FALSE
)

pbmc <- RunPCA(
    pbmc,
    npcs = 30,
    verbose = FALSE
)

pbmc <- RunUMAP(
    pbmc,
    dims = 1:30,
    verbose = FALSE
)


# --------------------------------------------------
# 3. Reorganize cell-type annotations
# --------------------------------------------------
annotation <- as.character(pbmc$Cluster_Annotation_Merged)

# Merge Tregs into CD4 T cells and MAIT cells into CD8 T cells
# according to the original analysis grouping.
annotation[annotation == "T_reg"] <- "CD4_T"
annotation[annotation == "MAIT"] <- "CD8_T"

pbmc$Cluster_Annotation_Merged2 <- annotation

# Exclude cells with undetermined annotations.
pbmc2 <- subset(
    pbmc,
    subset = Cluster_Annotation_Merged2 != "Undetermined"
)

# Define cell-type order.
celltype_levels <- c(
    "B_Naive", "B_Plasma", "B_SM",
    "CD4_T", "CD8_T", "NK",
    "C_Monocyte", "I_Monocyte", "NC_Monocyte",
    "cDCs", "pDCs", "Platelet",
    "Erythrocytes", "HSCs", "VD2p"
)

# Check for annotations missing from the defined levels.
other_types <- setdiff(
    unique(as.character(pbmc2$Cluster_Annotation_Merged2)),
    celltype_levels
)

if (length(other_types) > 0) {
    warning(
        "Additional cell types: ",
        paste(other_types, collapse = ", ")
    )
}

Idents(pbmc2) <- factor(
    pbmc2$Cluster_Annotation_Merged2,
    levels = c(celltype_levels, other_types)
)


# --------------------------------------------------
# 4. Visualize cell-type distribution using UMAP
# --------------------------------------------------
celltype_colors <- c(
    "B_Naive" = "#C04140",
    "B_Plasma" = "#ED7374",
    "B_SM" = "#CC85B1",
    "CD4_T" = "#7179A9",
    "CD8_T" = "#A5BCD7",
    "NK" = "#755A91",
    "C_Monocyte" = "#F2C5A7",
    "I_Monocyte" = "#F2DEC4",
    "NC_Monocyte" = "#F7E0D2",
    "cDCs" = "#FDFBE4",
    "pDCs" = "#F7E0D2",
    "Platelet" = "#B4B4B4",
    "Erythrocytes" = "#D6D6D6",
    "HSCs" = "#9CCDAA",
    "VD2p" = "#81B7A9"
)

p1 <- DimPlot(
    pbmc2,
    reduction = "umap"
) +
    scale_color_manual(values = celltype_colors) +
    theme(aspect.ratio = 1)

print(p1)


# --------------------------------------------------
# 5. Visualize T-cell marker expression
# --------------------------------------------------
p2 <- FeaturePlot(
    pbmc2,
    features = c("CD3E", "CD8A"),
    ncol = 1,
    keep.scale = "feature",
    pt.size = 1.5,
    cols = c("#EFECDC", "#C15858"),
    max.cutoff = 0.5
) +
    theme(aspect.ratio = 1)

print(p2)


# --------------------------------------------------
# 6. Identify DEGs between T1D and healthy controls
# --------------------------------------------------

# Set disease condition as the active identity.
Idents(pbmc2) <- "COND"

# CD4 T-cell differential expression.
cd4_obj <- subset(
    pbmc2,
    subset = Cluster_Annotation_Merged2 == "CD4_T"
)

deg_cd4 <- FindMarkers(
    cd4_obj,
    ident.1 = "T1D",
    ident.2 = "H",
    assay = "RNA",
    test.use = "wilcox",
    only.pos = FALSE
)

saveRDS(
    deg_cd4,
    "CD4T_disease_DEG_new.rds"
)


# CD8 T-cell differential expression.
cd8_obj <- subset(
    pbmc2,
    subset = Cluster_Annotation_Merged2 == "CD8_T"
)

deg_cd8 <- FindMarkers(
    cd8_obj,
    ident.1 = "T1D",
    ident.2 = "H",
    assay = "RNA",
    test.use = "wilcox",
    only.pos = FALSE
)

saveRDS(
    deg_cd8,
    "CD8T_disease_DEG_new.rds"
)

