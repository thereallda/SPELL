# monocle3 
# assume have processed seurat object with annotated celltypes
library(monocle3)
cds1 <- SeuratWrappers::as.cell_data_set(seu)
cds1 <- cluster_cells(cds1) # larger resolution represent more clusters
cds1@clusters$UMAP$partitions[cds1@clusters$UMAP$partitions == "2"] <- "1"

cds1 <- learn_graph(cds1, use_partition = F, verbose = FALSE)
plot_cells(cds1,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = F)
cds1 <- order_cells(cds1, root_cells = rownames(subset(colData(cds1), celltype == 'iPSC')))

plot_cells(cds1,
           color_cells_by = "pseudotime",
           group_cells_by = "celltype",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = F,
           trajectory_graph_color = "grey60")

tmp1 <- pseudotime(cds1)
write.csv(data.frame("monocle3_pseudotime"=tmp1), 
          file='results/trajectory/monocle3_pseudotime.csv')

seu <- AddMetaData(seu, tmp1, col.name = 'monocle3_pseudotime')
seu$monocle3_pseudotime_norm <- scale(log1p(seu$monocle3_pseudotime))
fp1 <- FeaturePlot(seu, "monocle3_pseudotime",
                   order=F, pt.size=.25) +
  scale_color_gradientn(colours = rev(c("#000004",'#4C2E7A',"#6F266D",'#BC5882','#D06E65','#D06E65','#EEEB89')))
ggsave("results/trajectory/Monocle3_pseudotime.png", fp1, width=6, height=5)
