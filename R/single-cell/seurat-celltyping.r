# lda's celltyping procedures 
library(Seurat)
library(tidyverse)
seu <- readRDS('seurat.rds')

# color 
# celltype
ct2col1 <- c("EC"="#FFCC48",
             "iPSC"="#66C2A5", 
             "MESD"="#F6883D",
             "MESC"="#1973AD"
             )
seu <- FindClusters(seu, resolution = 2)

dp2 <- DimPlot(seu, reduction="umap", 
               label=T, repel=T, ncol=4, order = T,
               group.by=c("seurat_clusters"), split.by = 'day')
ggsave("results/qc/umap_res2_splitday.pdf", dp2, width = 20, height = 8)

dp1 <- DimPlot(seu, reduction = "umap", 
               cols = c("grey","#A83C1E"), group.by = "genotype", split.by = 'day')
ggsave("results/qc/umap_genotype_day.png", dp1, width = 10, height = 4)

# cell markers ----
# to be refined
mrks1 <- c(
  'POU5F1','SOX2','NANOG','NCOA5','ESRP1','SPN','ZNF296', # iPSC
  'TBX6','MSGN1','FOXF2','WNT3', #Mesoderm
  'COL1A1','PDGFRB','TAGLN','SAMD3','BMP4','BMP5', #Mesenchymal
  'CDH5','PECAM1','DLL4','HEY2' #EC
)
# vis ----
## dotplot
dotp1 <- DotPlot(seu, features = mrks1, 
                 cols = c('#24bfca','#ee675b'),
                 group.by = "seurat_clusters") + coord_flip()
ggsave("results/celltype/cellmrks_dotplot_res2.png", dotp1, width=12, height=8, bg="white")

## featureplot
fp1 <- FeaturePlot(seu,
                   features = mrks1,
                   pt.size=.25, reduction = "umap", order = T, ncol = 3,
                   min.cutoff = "q10", max.cutoff = "q99") &
  scale_color_gradientn(colours = c("#cccccc","#fecb4c","#FEB500","#d34b26","#B41D6C"))
ggsave("results/celltype/cellmrks_umap_all_res2.pdf", fp1, width=15, height=30)
ggsave("results/celltype/cellmrks_umap_all_res2.png", fp1, width=15, height=30)
 
round(100*as.matrix(table(seu$seurat_clusters, seu$orig.ident))/as.vector(table(seu$orig.ident)),3)
DimPlot(seu, group.by = c("seurat_clusters","celltype"), label=T, repel=T)

# assign celltype ----
# MESD: Mesoderm; MESC: Mesenchymal; EC: endothelial cell
seu <- RenameIdents(seu, 
                                   "0" = "A", 
                                   "1" = "B", 
                                   "2" = "C"
)
seu$celltype <- Idents(seu)
seu$celltype <- factor(seu$celltype, levels=c('A','B','C'))
# save 
saveRDS(seu, file='results/int/seurat_celltype.rds')


DimPlot(seu, reduction = "umap", cols = ct2col1,
        group.by = "celltype", label=T, repel = T)
dp1 <- DimPlot(seu, reduction = "umap", cols = ct2col1,
               group.by = "celltype", label=T, repel = T)
ggsave("results/celltype/celltype_res2.png", dp1, width = 6, height = 5)
ggsave("results/celltype/celltype_res2.pdf", dp1, width = 6, height = 5)

dp2 <- DimPlot(seu, reduction = "umap", cols = ct2col1,
               group.by = "celltype", split.by = "day",label = F, repel = T)
ggsave("results/celltype/celltype_res2_byday.pdf", dp2, width = 16, height = 6)
ggsave("results/celltype/celltype_res2_byday.png", dp2, width = 10, height = 4)

dotp2 <- DotPlot(seu, features = mrks1, group.by = "celltype") + 
  theme(axis.text.x = element_text(angle=90)) + 
  scale_color_gradientn(colors = c("#cccccc","#fecb4c","#FEB500","#d34b26","#B41D6C")) +
  labs(x='',y='')
ggsave("results/celltype/cellmrks_celltype_dotplot_res2.png", dotp2, width=10, height=3.5, bg="white")

# heatmap
# average expression
gene_cell_exp <- AverageExpression(seu,
                                   features = mrks1,
                                   group.by = 'celltype',
                                   slot = 'data', assays = "RNA") 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

## heatmap ----
library(ComplexHeatmap)

marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
summary(marker_exp)

col_breaks <- c(-1,-0.5,0,0.5,1,1.5)
my_pal <- c("#053264","#569EC9","#FFFFFF","#F8B597","#C13438","#6C0321")
col_fun <- circlize::colorRamp2(col_breaks, my_pal)
pdf("results/celltype/cellmrks_heatmap.pdf", width=5, height=8)
hm1 <- Heatmap(marker_exp[mrks1,],
               name = "Z-score", 
               cluster_columns = F, cluster_rows = F,
               column_names_rot = 90,
               show_column_names = T, show_row_names = T,
               col = col_fun, use_raster = F
)
draw(hm1)
dev.off()

# basic cell composition in each group or sample ----
seu@meta.data[1:3,]
# by group
df1 <- data.frame(table(seu$orig.ident, seu$celltype))
df1$day <- str_split(df1$Var1, '_', simplify = T)[,1]
df1$genotype <- str_split(df1$Var1, '_', simplify = T)[,2]

# number of celltypes
bp1 <- ggplot(df1, aes(Var2, Freq)) +
  geom_bar(aes(fill=Var2),stat="identity", width=0.6) +
  geom_text(aes(label=Freq),hjust = -1, 
            data=df1 %>% group_by(Var2) %>% summarise(Freq=sum(Freq)) %>% filter(Freq<5000)) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.position = 'none') +
  scale_fill_manual(values = ct2col1) +
  labs(x='', y='Number of Cells')
ggsave('results/celltype/CellNum_celltype.png', bp1, width=6, height=3.5)

# numbers of celltypes by day
bp3 <- ggplot(df1, aes(Var2, Freq)) +
  geom_bar(aes(fill=genotype),stat="identity", width=0.6, 
           position = position_dodge()) +
  facet_wrap(~day,nrow=1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        strip.background = element_blank()) +
  scale_fill_manual(values=c("grey","#A83C1E")) +  
  labs(x='', y='Number of Cells',fill='')

# Proportion of celltypes by day 
df2 <- df1 %>% 
  group_by(Var1) %>% 
  mutate(ncell_day=sum(Freq),
         prop_day=round(100*Freq/ncell_day,3)) %>% 
  arrange(day) 

bp4 <- ggplot(df2, aes(Var2, prop_day)) +
  geom_bar(aes(fill=genotype),stat="identity", width=0.6, 
           position = position_dodge()) +
  facet_wrap(~day,nrow=1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        strip.background = element_blank()) +
  scale_fill_manual(values=c("grey","#A83C1E")) +  
  labs(x='', y='Percentage of Cells (%)',fill='')
bps1 <- bp3/bp4
ggsave('results/celltype/CellProp_day_celltype.png', bps1, width=10, height=8)
