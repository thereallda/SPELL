## peak annotation
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

peaks <- read.table('test.bed')
head(peaks)
# V1      V2      V3    V4       V5
# 1 chr1  629803  630074 peak1 61.28789
# 2 chr1  633893  634164 peak2 64.13491
# 3 chr1  778616  778859 peak3  4.77470
# 4 chr1  958827  958973 peak4  4.99200
# 5 chr1 1001089 1001296 peak5  7.56939
# 6 chr1 1006450 1006883 peak6 19.31422

# Convert to GRanges
gr1 <- GRanges(
  seqnames = Rle(peaks$V1),
  ranges = IRanges(start = peaks$V2, end = peaks$V3),
  score = peaks$V5
)

# Annotate the peaks using the RefSeq-based TxDb object
peakAnno1 <- annotatePeak(gr1, 
                          TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, 
                          annoDb = "org.Hs.eg.db",  # Use org.Hs.eg.db to add gene symbols and names
                          verbose = FALSE)

pdf('results/peaks_distribution.pdf', width=6, height=4)
par(oma=c(0,0,2,0))
plotAnnoPie(peakAnno1, main='ChIP-seq peaks')
dev.off()

# save table
peaks_anno1 <- as.data.frame(peakAnno1)
write.csv(peaks_anno1, file='peakAnnotation.csv', row.names = F)
