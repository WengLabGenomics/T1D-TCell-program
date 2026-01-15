library(dplyr)
library(Seurat)
library(ggplot2)

# -----------------------------
# Load Seurat object
# -----------------------------
obj <- readRDS("T1D_Seurat_Object_Final.rds")

# -----------------------------
# Dimensionality reduction 
# -----------------------------
pbmc <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

# Exclude cells annotated as "Undetermined" based on the original study’s cell-type labels
pbmc2 <- subset(pbmc, subset = Cluster_Annotation_Merged != "Undetermined")

# -----------------------------
# UMAP visualization (manual colors)
# -----------------------------
celltype_colors <- c(
  "B_Naive"      = "#8B0000",
  "B_Plasma"     = "#FF69B4",
  "B_SM"         = "#CD5C5C",
  "C_Monocyte"   = "#FFDAB9",
  "CD4_T"        = "#0000EE",
  "CD8_T"        = "#00F5FF",
  "NK"           = "#AB82FF",
  "I_Monocyte"   = "#CD6600",
  "NC_Monocyte"  = "#8B4500",
  "cDCs"         = "#E9967A",
  "pDCs"         = "#FF8C00",
  "Platelet"     = "#696969",
  "Erythrocytes" = "grey"
)

# If you want to save the plot, uncomment the pdf/dev.off lines
# pdf("UMAP.pdf", width = 7, height = 6)
p_umap <- DimPlot(pbmc2, reduction = "umap", group.by = "Cluster_Annotation_Merged") +
  scale_color_manual(values = celltype_colors) +
  theme(aspect.ratio = 1)
print(p_umap)
# dev.off()

# -----------------------------
# DotPlot of canonical marker genes
# -----------------------------
cd_genes <- c(
  "CD3E", "CD40LG", "CD8A", "TRDV2", "NCAM1",
  "CD79A", "MS4A1", "XBP1",
  "CD14", "FCGR3A",
  "FCER1A", "LILRA4",
  "CD34", "PPBP", "HBB"
)

# pdf("dotplot-marker.pdf", width = 10, height = 6)
p_dot <- DotPlot(
  pbmc2,
  features = cd_genes,
  group.by = "Cluster_Annotation_Merged",
  scale = TRUE,
  scale.by = "size"
) +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(
      angle = 90,     # Rotate gene labels for readability
      hjust = 1,
      vjust = 0.5
    )
  )
print(p_dot)
# dev.off()

# -----------------------------
# Differential expression (T1D vs healthy) within CD4_T and CD8_T
# -----------------------------
# Set identity class to condition label
Idents(pbmc2) <- "COND"

# CD4 T cells
cd4_obj <- subset(pbmc2, subset = Cluster_Annotation_Merged %in% "CD4_T")
cd4_markers <- FindMarkers(cd4_obj, ident.1 = "T1D", ident.2 = "H")
saveRDS(cd4_markers, "CD4T_DEG.rds")

# CD8 T cells
cd8_obj <- subset(pbmc2, subset = Cluster_Annotation_Merged %in% "CD8_T")
cd8_markers <- FindMarkers(cd8_obj, ident.1 = "T1D", ident.2 = "H")
saveRDS(cd8_markers, "CD8T_DEG.rds")