library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(tidyverse)
##--Custom Function--##
barplot_pubr <- function(enrich.obj, showCategory=10){
  eobj <- enrich.obj[,c('Description','Count','p.adjust')]
  eggtmp <- eobj %>% 
    head(.,n=showCategory) %>% 
    dplyr::mutate(p.adjust = -log10(p.adjust)) %>% 
    reshape2::melt()
  bp <- ggplot(eggtmp,aes(x=factor(Description, levels=rev(unique(Description))))) +
    geom_bar(aes(y=value, fill=variable), stat='identity', 
             position=position_dodge(width=0.9), width=0.8) +
    coord_flip() +
    theme_minimal() +
    geom_hline(yintercept = (-log10(0.05)),lty=4,col="grey",lwd=0.6) +
    theme(axis.text.y = element_text(size=15, color='black'),
          axis.text.x = element_text(color='black'),
          axis.line = element_line(),
          axis.ticks = element_line(color='black'),
          legend.position = 'top',
          panel.grid = element_blank()) +
    scale_fill_manual(values = c("#8696A1","#CDDEE9"), labels = c('Gene Counts', '-log10(p.adjust)')) +
    labs(x='',y='',fill='')
  return(bp)
}
##----##
ego1 <- enrichGO(gene = de$GeneID,OrgDb = org.Mm.eg.db,keyType = 'ENSEMBL',ont = 'BP')
bp1 <- barplot_pubr(ego1)
ggsave('results/GOBP_v1.pdf',bp1,width=8,height=3)
# simplify redudant GO terms, not working if `ont="ALL"`
ego1.sim <- simplify(ego1, cutoff=0.6)
bp2 <- barplot_pubr(ego1.sim)


# compare enrichment
expr2num <- function(expr.char){
  num <- as.numeric(str_split(expr.char,'/',simplify = T)[,1])/as.numeric(str_split(expr.char,'/',simplify = T)[,2])
}
compare_dotplot <- function(enrich.compare.obj, showCategory=10){
  eobj <- as.data.frame(enrich.compare.obj)
  eobj %>% 
    dplyr::group_by(Cluster) %>% 
    slice_min(order_by = p.adjust, n=showCategory) %>% 
    ungroup() %>% 
    dplyr::select(Description) %>% 
    dplyr::inner_join(eobj, y=., by='Description') %>%
    ggplot(aes(x=Cluster, 
               y=factor(Description, levels=rev(unique(Description))))) +
    geom_point(aes(size=expr2num(GeneRatio), color=p.adjust)) + 
    geom_text(aes(label=paste0(round(expr2num(GeneRatio),2)*100, '%')),nudge_x = 0.15) + 
    scale_color_gradient(low="#f0c27b",high="#4b1248") + 
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.ticks = element_line(color='black'),
          axis.text = element_text(color='black', size=12)) +
    labs(x='',y='',size='GeneRatio')
}

##----##
egos3 <- compareCluster(list(cluster1=de1$GeneID,
                             cluster2=de2$GeneID,
                             cluster3=de3$GeneID),
                        fun='enrichGO', OrgDb='org.Mm.eg.db',keyType='ENSEMBL',ont='BP')
egos3 <- setReadable(egos3, org.Mm.eg.db)
dp3 <- compare_dotplot(egos3, showCategory=10)
ggsave('results/Compare_GOBP_v1.pdf', dp3, width=7, height=3)
