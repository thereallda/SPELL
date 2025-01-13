library(edgeR)
library(tidyverse)
# make directory for results
if (!dir.exists('results/qc')) dir.create('results/qc', recursive = T)
if (!dir.exists('results/de')) dir.create('results/de', recursive = T)
# read in meta information (required)
metadata <- read.csv("metadata.csv", header=TRUE)

# step 0. Merging counts from featureCounts ----
counts_df <- metadata$path %>%
  map_dfc(., ~ read.table(.x, sep = '\t', row.names = 1, skip = 2,
                          colClasses = c('character',rep('NULL',5),'integer')))

colnames(counts_df) <- metadata$id
# if exists, remove version number in ENSEMBL gene id 
#rownames(counts_df) <- gsub('\\..*', '', rownames(counts_df))
write.table(counts_df, file='data/Counts.csv', sep=',')

# step 1. filter lowly expressed genes ----
counts_keep <- counts_df[rowSums(counts_df>0) > 0.75*ncol(counts_df), ]
# step 2. DE and output de results with edgeR ----
samples_group <- metadata$condition

degs <- DGEList(counts=counts_keep, group=samples_group)
degs <- calcNormFactors(degs)
counts_norm <- cpm(degs)
write.table(counts_norm, file='data/Counts_norm.csv', sep=',')

bcv = 0.1 
# set up contrast
# for `exactTest`: if the pair is c("A","B") then the comparison is B - A,
contr.ls <- list(
    c('WT', 'MutA'),
    c('WT', 'MutB'),
    c('WT', 'MutC')
)
# extract DE results by contrasts
res.ls <- lapply(contr.ls, function(contr) {
    et <- exactTest(degs, dispersion=bcv^2, pair=contr)  
    res <- topTags(et, n=Inf)
})
names(res.ls) <- paste0('Res_', unlist(lapply(contr.ls, paste, collapse='_')))
# save unfiltered DE results
lapply(seq_along(res.ls), function(x, nm, i) {
      write.csv(x = x[[i]], file = paste0('results/de/', nm[[i]], '_DE.csv'))
}, x=res.ls, nm=names(res.ls))

# significant DEGs are considered as genes that pass cutoff:
# |fold change| >= 2 and p.adj < 0.05
res.sig.ls <- lapply(seq_along(res.ls), function(x, nm, i) {
  sig <- subset(x[[i]]$table, FDR < 0.05 & abs(logFC) >= 1)
  # save significant DEGs tables
  write.csv(x = sig, file = paste0('results/de/', nm[[i]], '_DE_sig.csv'), row.names=FALSE)
  return(sig)
}, x=res.ls, nm=names(res.ls))
# save DE results as .RData
save(de, res.ls, res.sig.ls, file = 'data/DE.RData')