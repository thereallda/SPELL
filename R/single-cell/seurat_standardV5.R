install.packages('Seurat')
library(Seurat)

seu.counts <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")
seu <- CreateSeuratObject(counts = seu.counts)
# quality control ----
# adjust criteria based on data
seu <- subset(seu, subset = nCount_RNA > 150 & nCount_RNA < 3500 &
                nFeature_RNA > 150 & nFeature_RNA < 1500)
seu <- NormalizeData(object = seu)
seu <- FindVariableFeatures(object = seu, nfeatures = nrow(seu))
seu <- ScaleData(object = seu, features=rownames(seu)) # use all genes
seu <- RunPCA(object = seu)
# check appropriate numbers of PCs, such as 1:20
ElbowPlot(seu)
seu <- FindNeighbors(object = seu, dims = 1:20)
# resolution can be tuned based on the clustering effect
seu <- FindClusters(object = seu, resolution = 0.8)
seu <- RunUMAP(object = seu)
seu <- RunTSNE(object = seu, dims = 1:20)
DimPlot(object = seu, reduction = "tsne")
# find markers for every cluster compared to all remaining cells, report only the positive ones
seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# cell typing ----
# set resolution = 0.8
Idents(seu) <- "RNA_snn_res.0.8"
# annotate cell type by cluster
seu <- RenameIdents(seu,
                         `0` = "Muscle Cell",
                         `1` = "Muscle Cell",
                         `2` = "Epithelial Cell",
                         `3` = "Neuron",
                         `4` = "Unknown",
                         `5` = "Glia",
                         `6` = "Tracheal Cell",
                         `7` = "Neuron",
                         `8` = "Epithelial Cell",
                         `9` = "Fat Body",
                         `10` = "Glia",
                         `11` = "Unknown",
                         `12` = "Unknown",
                         `13` = "Fat Body",
                         `14` = "Unknown",
                         `15` = "Germline Cell",
                         `16` = "Unknown"
  
)
# plot by genes
FeaturePlot(object = seu, features = 'Gsdmd')
