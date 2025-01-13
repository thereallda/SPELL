# integration with Seurat v5
library(Seurat)
library(tidyverse)
library(patchwork)
load("seurat.RData")
count_mat <- seu.combined@assays$RNA$counts
meta <- seu.combined@meta.data

seu <- CreateSeuratObject(counts = count_mat, meta.data = meta)
unique(seu$orig.ident)
# [1] "O1" "O2" "Y1" "Y2"

# split different batch into different layers
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
# Note that since the data is split into layers, normalization and variable 
# feature identification is performed for each batch independently
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50)
seu <- FindNeighbors(seu, dims = 1:30, reduction = "pca")
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca")
DimPlot(seu, reduction = "umap", group.by = c("orig.ident"))
ggsave("results/qc/UMAP_Unintegrated.png", width=12, height=6)

# Perform streamlined (one-line) integrative analysis
# CCA
seu <- IntegrateLayers(seu, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = T
)
seu <- FindNeighbors(seu, reduction = "integrated.cca", dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- FindNeighbors(seu, reduction = "integrated.cca", dims = 1:30)
seu <- RunUMAP(seu, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap")
p1 <- DimPlot(seu, reduction = "umap", group.by = c("orig.ident"))
ggsave("results/qc/UMAP_CCAintegrated.png", p1, width=12, height=6)

# Harmony
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = FALSE
)
seu <- FindNeighbors(seu, reduction = "integrated.harmony", dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- FindNeighbors(seu, reduction = "integrated.harmony", dims = 1:30)
seu <- RunUMAP(seu, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap")
p1 <- DimPlot(seu, reduction = "umap", group.by = c("orig.ident"))
ggsave("results/qc/UMAP_Harmonyintegrated.png", p1, width=12, height=6)

# FastMNN
seu <- IntegrateLayers(
  object = seu, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
seu <- FindNeighbors(seu, reduction = "integrated.mnn", dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- FindNeighbors(seu, reduction = "integrated.mnn", dims = 1:30)
seu <- RunUMAP(seu, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap")
p1 <- DimPlot(seu, reduction = "umap", group.by = c("orig.ident"))
ggsave("results/qc/UMAP_MNNintegrated.png", p1, width=12, height=6)
# for final integratio, rejoin data 
seu <- JoinLayers(seu)
save(seu, file="data/seu_integrated.RData")
