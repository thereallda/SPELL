# Basic DESeq2
library(DESeq2)
counts <- read.csv('Counts_filtered.csv', row.names = 1)
counts_keep <- counts[rowSums(counts > 1) > 0.7*ncol(counts), ]

colData1 <- data.frame(row.names = colnames(counts_keep),
                   condition = rep(c('A','B'), each=8))
colData1
#     condition
# S1          A
# S2          A
# S3          A
# S4          A
# S5          A
# S6          A
# S7          A
# S8          A
# S9          B
# S10         B
# S11         B
# S12         B
# S13         B
# S14         B
# S15         B
# S16         B

dds <- DESeqDataSetFromMatrix(countData = counts_keep,
                              colData = colData1,
                              design = ~0+condition)
dds <- DESeq(dds)

deres1 <- results(dds, contrast = c('condition','A','B'), tidy=TRUE)
head(deres1)
#                   row   baseMean log2FoldChange      lfcSE        stat       pvalue
# 1 ENSMUSG00000025902   16.08081     0.31224893 0.27098405   1.1522779 2.492069e-01
# 2 ENSMUSG00000098104   11.50252    -0.51520680 0.31064337  -1.6585153 9.721349e-02
# 3 ENSMUSG00000103922  428.29390    -1.46714681 0.11895438 -12.3336935 5.964858e-35
# 4 ENSMUSG00000033845  601.32767     0.08623021 0.09793307   0.8805014 3.785877e-01
# 5 ENSMUSG00000025903 1192.30860    -0.26643769 0.08696370  -3.0637805 2.185592e-03
# 6 ENSMUSG00000033813  176.99666     0.53226849 0.13242355   4.0194399 5.833666e-05
#           padj
# 1 3.176612e-01
# 2 1.386648e-01
# 3 1.207884e-33
# 4 4.533104e-01
# 5 4.347064e-03
# 6 1.437774e-04

# if have multiple contrast group
colData2 <- data.frame(
  row.names = colnames(counts_keep),
  condition = rep(c('A','B','C','D'), each=4)
)

contrast.df <- data.frame(
  Group1 = c('A','A','A','B','B'),
  Group2 = c('B','C','D','C','D')
)

dds <- DESeqDataSetFromMatrix(countData = counts_keep,
                              colData = colData2,
                              design = ~0+condition)
dds <- DESeq(dds)

res.ls <- lapply(seq_len(nrow(contrast.df)), function(i) {
  restmp <- results(dds, contrast = c('condition', contrast.df[i,1],contrast.df[i,2]), tidy=TRUE) %>% 
    dplyr::rename(GeneID = row)
})
names(res.ls) <- paste0('Res_', apply(contrast.df, 1, paste, collapse='_'))

# save unfiltered DE results
lapply(seq_along(res.ls), function(x, nm, i) {
  write.csv(x = x[[i]], file = paste0(outdir, nm[[i]], '_DE.csv'), row.names=FALSE)
}, x=res.ls, nm=names(res.ls))
