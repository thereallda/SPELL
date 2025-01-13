# seurat find markers and enrichment analysis
library(Seurat)
library(tidyverse)
library(gprofiler2)

# barplot dot
lollipopPlot <- function(enrich.obj, 
                         term_col="Description", 
                         term_size_col="Count", 
                         p_col="p.adjust",
                         showCategory=10, title=NULL, str_width=35, 
                         log=TRUE,
                         palette='#689C9A'){
  eobj <- enrich.obj[, c(term_col, term_size_col, p_col)]
  
  eggtmp <- head(eobj, n=showCategory)
  
  if (log) {
    eggtmp[, p_col] <- -log10(eggtmp[, p_col])
  }
  
  eggtmp[, term_col] <- factor(eggtmp[, term_col], levels=rev(unique(eggtmp[, term_col])))
  
  bp <- ggplot(eggtmp, aes(.data[[p_col]], .data[[term_col]])) +
    geom_bar(stat='identity', width=0.05, fill=palette) +
    geom_point(aes(size=.data[[term_size_col]]), color=palette) +
    theme_classic() +
    # geom_vline(xintercept = (-log10(0.05)),lty=4,col="grey",lwd=0.6) +
    theme(axis.text = element_text(color='black'),
          legend.position = 'right') +
    scale_y_discrete(labels = function(x) str_wrap(x, width = str_width)) +
    scale_x_continuous(expand = expansion(mult = c(0, .1))) +
    # scale_fill_manual(values = c("#2d74b9","#d28b47"), labels = c('Gene Counts', '-Log10(p.adjust)')) +
    labs(x='-log10(FDR)',y='',size='Enriched gene set size', title=title)
  return(bp)
}

# DEG for single cell data
seu$ct_group <- paste(seu$celltype, seu$group, sep='.')
table(seu$ct_group)

# DEGs
# celltype
ct1 <- c('CM','EC','EndoC')
deg_ls1 <- list()
for (i in ct1) {
  ct.ctri <- paste0(i, '.Ctrl')
  ct.koi <- paste0(i, '.Mut')
  degi <- FindMarkers(seu_integrated, 
                      ident.1 = ct.koi,
                      ident.2 = ct.ctri,
                      group.by = 'ct_group',
                      min.pct = 0.1)
  degi$GeneID <- rownames(degi)
  degi_sig <- subset(degi, p_val_adj < 0.05)
  deg_ls1[[paste0(ct.koi, '_vs_', ct.ctri)]] <- degi
  # save table
  write.csv(degi, file=paste0('results/de/DEGs_',i,'_MutvsCtrl.csv'))
  write.csv(degi_sig, file=paste0('results/de/DEGs_',i,'_LMutvsCtrl_sig.csv'))
}

# number of DEG per group

# enrichment ----
## single celltype ----
df1 <- subset(deg_ls1[[1]], p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5))
ego1 <- gost(df1$GeneID, exclude_iea = T, evcodes = T)
gostplot(ego1)

# top enrichment term
edf1 <- ego1$result %>% 
  arrange(p_value) %>% 
  filter(source %in% c("GO:BP","KEGG","REAC")) 
lp1 <- lollipopPlot(edf1, p_col = 'p_value', 
                    term_col = 'term_name', term_size_col = 'intersection_size',
                    title='DEGs (Mutant vs Control)')
## all celltypes ----
# 