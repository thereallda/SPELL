# gprofiler enrichment
library(gprofiler2)
# custom function
barplot_pubr <- function(enrich.obj, 
                         term_col="Description", 
                         term_size_col="Count", 
                         p_col="p.adjust",
                         showCategory=10, title=NULL, str_width=35, 
                         log=TRUE){
  eobj <- enrich.obj[, c(term_col, term_size_col, p_col)]
  
  eggtmp <- head(eobj, n=showCategory)
  
  if (log) {
    eggtmp[, p_col] <- -log10(eggtmp[, p_col])
  }
  
  eggtmp <- reshape2::melt(eggtmp)
  eggtmp[, term_col] <- factor(eggtmp[, term_col], levels=rev(unique(eggtmp[, term_col])))
  
  bp <- ggplot(eggtmp, aes_string(term_col)) +
    geom_bar(aes(y=value, fill=variable), stat='identity', 
             position=position_dodge(width=0.9), width=0.8) +
    coord_flip() +
    theme_minimal() +
    geom_hline(yintercept = (-log10(0.05)),lty=4,col="grey",lwd=0.6) +
    theme(axis.text = element_text(color='black'),
          axis.line = element_line(),
          axis.ticks = element_line(color='black'),
          legend.position = 'right',
          panel.grid = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = str_width)) +
    scale_fill_manual(values = c("#2d74b9","#d28b47"), labels = c('Gene Counts', '-Log10(p.adjust)')) +
    labs(x='',y='',fill='', title=title)
  return(bp)
}

goM1 <- gost(de$geneID,
     organism = 'hsapiens', exclude_iea = TRUE, 
     correction_method = 'fdr', evcodes = TRUE)
gostplot(goM1)

# top enrichment term
edf1 <- goM1$result %>% 
  arrange(p_value) %>% 
  filter(source %in% c("GO:BP","KEGG","REAC")) 

barplot_pubr(edf1$result, 
            p_col = 'p_value', term_col = 'term_name', term_size_col = 'intersection_size')
