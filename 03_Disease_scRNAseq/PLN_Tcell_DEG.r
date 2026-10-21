
library(Seurat)
library(dplyr)
library(scDblFinder)
library(BiocParallel)
library(harmony)


# --------------------------------------------------
# 1. Load and merge PLN CITE-seq samples
# --------------------------------------------------
samples <- c(
    "HPAP012", "HPAP020", "HPAP021", "HPAP023",
    "HPAP026", "HPAP027", "HPAP032", "HPAP048",
    "HPAP095", "HPAP098", "HPAP099", "HPAP102",
    "HPAP110"
)

for (i in samples) {

    data <- Read10X(
        paste0(
            "/home/lifesci/shiyu0724/workspace/eGene-Knock/PANC/CITE/hpapdata/",
            i, "_count/outs/filtered_feature_bc_matrix"
        )
    )

    obj <- CreateSeuratObject(
        counts = data[["Gene Expression"]],
        project = i
    )

    obj[["ADT"]] <- CreateAssayObject(
        counts = data[["Antibody Capture"]]
    )

    obj$sample <- i

    new_barcodes <- paste0(i, "_", colnames(obj))
    obj <- RenameCells(obj, new.names = new_barcodes)

    if (i == "HPAP012") {
        obj_all <- obj
    } else {
        obj_all <- merge(obj_all, y = obj)
    }
}

saveRDS(obj_all, "raw.rds")


# --------------------------------------------------
# 2. Detect doublets and perform quality control
# --------------------------------------------------
DefaultAssay(obj_all) <- "RNA"
obj_all[["RNA"]] <- JoinLayers(obj_all[["RNA"]])

sce <- as.SingleCellExperiment(obj_all, assay = "RNA")

sce <- scDblFinder(
    sce,
    samples = "sample",
    clusters = FALSE,
    BPPARAM = MulticoreParam(10)
)

obj <- obj_all

obj$scDblFinder.score <- colData(sce)$scDblFinder.score
obj$scDblFinder.class <- colData(sce)$scDblFinder.class

# Retain singlets and high-quality cells.
obj <- subset(
    obj,
    subset = scDblFinder.class == "singlet"
)

obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
)

obj <- subset(
    obj,
    subset = nFeature_RNA > 500 &
             nFeature_RNA < 8500 &
             percent.mt < 15
)


# --------------------------------------------------
# 3. Normalize RNA and integrate samples using Harmony
# --------------------------------------------------
DefaultAssay(obj) <- "RNA"

# Split RNA layers by sample for Seurat v5 integration.
obj[["RNA"]] <- split(
    obj[["RNA"]],
    f = obj$sample
)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

obj <- IntegrateLayers(
    object = obj,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    assay = "RNA",
    verbose = FALSE
)


# --------------------------------------------------
# 4. Normalize ADT and integrate samples using Harmony
# --------------------------------------------------
DefaultAssay(obj) <- "ADT"

# Use all ADT features for dimensionality reduction.
VariableFeatures(obj) <- rownames(obj[["ADT"]])

obj <- NormalizeData(
    obj,
    normalization.method = "CLR",
    margin = 2
) %>%
    ScaleData() %>%
    RunPCA(reduction.name = "apca")

obj <- RunHarmony(
    object = obj,
    group.by.vars = "sample",
    reduction = "apca",
    reduction.save = "aharmony"
)


# --------------------------------------------------
# 5. Perform WNN analysis and clustering
# --------------------------------------------------

# Combine integrated RNA and ADT embeddings.
obj <- FindMultiModalNeighbors(
    obj,
    reduction.list = list("harmony", "aharmony"),
    dims.list = list(1:30, 1:20),
    modality.weight.name = "RNA.weight"
)

obj <- RunUMAP(
    obj,
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_"
)

obj <- FindClusters(
    obj,
    graph.name = "wsnn",
    algorithm = 3,
    resolution = 0.5,
    verbose = FALSE
)

saveRDS(obj, "wnn.rds")


# --------------------------------------------------
# 6. Visualize cell-type markers
# --------------------------------------------------

# ADT marker expression.
DefaultAssay(obj) <- "ADT"

FeaturePlot(
    obj,
    reduction = "wnn.umap",
    features = c(
        "Hu.CD3-UCHT1",
        "Hu.TCR.AB",
        "Hu.CD4-RPA.T4",
        "Hu.CD8",
        "Hu.CD2",
        "Hu.CD56",
        "Hu.CD19",
        "Hu.CD14-M5E2",
        "Hu.CD11c"
    )
)

# RNA marker expression.
DefaultAssay(obj) <- "RNA"

FeaturePlot(
    obj,
    reduction = "wnn.umap",
    features = c(
        "CD3D", "CD4", "CD8A",
        "MS4A1", "MZB1", "GNLY",
        "CD14", "FCGR3A", "FCER1A",
        "LYZ", "CLEC4C", "PPBP"
    )
)


# --------------------------------------------------
# 7. Annotate cell types and disease status
# --------------------------------------------------

# Assign cell-type labels based on RNA and ADT markers.
new.cluster.ids <- c(
    "B", "B", "B", "CD4T", "CD4T", "B", "CD4T",
    "CD4T", "CD8T", "CD4T", "CD8T", "CD8T", "CD8T", "B",
    "B", "B", "B", "NK", "NK", "CD8T", "Monocyte",
    "B", "Plasmablast", "B", "B", "CD4T", "CD4T",
    "B", "B", "B", "B", "B", "B", "B"
)

# Confirm that the number of labels matches the clusters.
stopifnot(
    length(new.cluster.ids) == length(levels(obj))
)

names(new.cluster.ids) <- levels(obj)

obj <- RenameIdents(obj, new.cluster.ids)
obj$cell_type <- Idents(obj)

# Define T1D and healthy control samples.
t1d <- c(
    "HPAP020", "HPAP021", "HPAP023",
    "HPAP032", "HPAP048", "HPAP098", "HPAP102"
)

hc <- c(
    "HPAP012", "HPAP026", "HPAP027",
    "HPAP095", "HPAP099", "HPAP110"
)

obj$Disease <- NA_character_

obj$Disease[obj$sample %in% t1d] <- "T1D"
obj$Disease[obj$sample %in% hc] <- "HC"


# --------------------------------------------------
# 8. Identify DEGs in CD8 T cells
# --------------------------------------------------
cd8 <- subset(
    obj,
    subset = cell_type == "CD8T"
)

Idents(cd8) <- "Disease"
DefaultAssay(cd8) <- "RNA"

cd8[["RNA"]] <- JoinLayers(cd8[["RNA"]])

deg_cd8 <- FindMarkers(
    cd8,
    ident.1 = "T1D",
    ident.2 = "HC",
    logfc.threshold = 0,
    test.use = "wilcox"
)

write.csv(
    deg_cd8,
    "CD8T_DEG_LD.csv"
)


# --------------------------------------------------
# 9. Identify DEGs in CD4 T cells
# --------------------------------------------------
cd4 <- subset(
    obj,
    subset = cell_type == "CD4T"
)

Idents(cd4) <- "Disease"
DefaultAssay(cd4) <- "RNA"

cd4[["RNA"]] <- JoinLayers(cd4[["RNA"]])

deg_cd4 <- FindMarkers(
    cd4,
    ident.1 = "T1D",
    ident.2 = "HC",
    logfc.threshold = 0,
    test.use = "wilcox"
)

write.csv(
    deg_cd4,
    "CD4T_DEG_LD.csv"
)
