# Basic edgeR
library(edgeR)
counts <- read.csv('Counts_filtered.csv', row.names = 1)
counts_keep <- counts[rowSums(counts > 1) > 0.7*ncol(counts), ]

meta <- data.frame(
  row.names = colnames(counts_keep),
  condition = rep(c('A','B','C','D'), each=4)
)
contrast.df <- data.frame(
  Group1 = c('A','A','A','B','B'),
  Group2 = c('B','C','D','C','D')
)

degs <- edgeR::DGEList(counts, group = meta$condition)
degs <- edgeR::calcNormFactors(degs, method = "TMM") # Default: perform TMM normalization
design.mat <- model.matrix(~0+condition, data = meta)
degs <- edgeR::estimateDisp(degs, design = design.mat)
fit.glm <- edgeR::glmFit(degs, design = design.mat)

# construct contrast
contrast.vec <- apply(contrast.df, 1, function(x) { paste(paste0("condition",x), collapse="-") })
de.ls2 <- lapply(1:length(contrast.vec), function(i) {
  contrast.mat <- limma::makeContrasts(contrasts = contrast.vec[i], levels=design.mat)
  # LRT 
  lrt.glm <- edgeR::glmLRT(fit.glm, contrast = contrast.mat)
  
  # extract DE results
  res1 <- edgeR::topTags(lrt.glm, n = Inf, adjust.method = "BH")
  res.tab <- res1$table %>% tibble::rownames_to_column(var = "GeneID")
  
  # filter significant DEGs
  res.sig.tab <- res.tab[abs(res.tab$logFC) >= logfc.cutoff & res.tab$FDR < p.cutoff,]
  return(list(res.tab = res.tab, res.sig.tab = res.sig.tab))
})
# unlist 
res.ls <- lapply(de.ls2, function(x) x$res.tab)
res.sig.ls <- lapply(de.ls2, function(x) x$res.sig.tab)
names(res.ls) <- names(res.sig.ls) <- gsub("-","_",gsub("condition","",contrast.vec))

