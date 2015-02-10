# TOOL cna-correlate.R: "Plot correlations of called copy number data" (Calculates correlations between identified regions in called copy number data and plots them. This is recommended to be performed after running the Identify common regions from called copy number data tool.)
# INPUT regions.tsv: regions.tsv TYPE GENERIC
# OUTPUT correlations.pdf: correlations.pdf
# OUTPUT correlations-scaled.pdf: correlations-scaled.pdf
# OUTPUT correlations.tsv: correlations.tsv
# PARAMETER genome.build: "Genome build" TYPE [GRCh37: GRCh37, NCBI36: NCBI36, NCBI35: NCBI35, NCBI34: NCBI34] DEFAULT GRCh37 (The genome build to use. GRCh37 = hg19, NCBI36 = hg18, NCBI35 = hg17, NCBI34 = hg16.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-25

source(file.path(chipster.common.path, 'library-Chipster.R'))
library(RColorBrewer)

dat <- readData('regions.tsv')

calls <- as.matrix(dat[,grep("^flag\\.", names(dat))])
calls[calls ==  2] <-  1
calls[calls == -2] <- -1

m <- cor(t(calls), method='kendall')

chrom <- chromosomeToInteger(dat$chromosome)
uni.chrom <- unique(chrom)

# bitmap('correlations.png', width=297, height=210, units="mm", res=150)
pdf('correlations.pdf', width=11.69, height=8.27)
image(x=1:(nrow(m)+1), y=1:(nrow(m)+1), z=m, zlim=c(-1,1), col=brewer.pal(10, 'RdBu'), xlab=NA, ylab=NA, useRaster=TRUE, main='Correlation Matrix', xaxt='n', yaxt='n')
chr.change <- chrom != c(chrom[-1], 0)
chr.change <- which(chr.change) + 1
abline(h=chr.change[-length(chr.change)], v=chr.change[-length(chr.change)], lty='dashed')
if (nrow(m) < 100 && 'cytoband' %in% colnames(dat)) {
  axis(side=1, at=1:nrow(m)+.5, labels=dat$cytoband, las=2, line=-.5, tick=FALSE, cex.axis=.4)
  axis(side=2, at=1:nrow(m)+.5, labels=dat$cytoband, las=2, line=-.5, tick=FALSE, cex.axis=.1)
  axis(side=4, at=1:nrow(m)+.5, labels=dat$cytoband, las=2, line=-.5, tick=FALSE, cex.axis=.4)
} else {
  ax <- (chr.change + c(0, chr.change[-length(chr.change)]))/2
  axis(side=1, at=ax, labels=uni.chrom, lwd=.5, las=1, cex.axis=1)
  axis(side=2, at=ax, labels=uni.chrom, lwd=.5, las=1, cex.axis=.5)
  axis(side=4, at=ax, labels=uni.chrom, lwd=.5, las=1, cex.axis=.5)
}
dev.off()

# bitmap('correlations-scaled.png', width=297, height=210, units="mm", res=150)
pdf('correlations-scaled.pdf', width=11.69, height=8.27)
pos <- dat$start
pos2 <- dat$end
chrom.lengths <- CGHbase:::.getChromosomeLengths(genome.build)[as.character(uni.chrom)]
chrom.ends <- integer()
cumul <- 0
for (j in uni.chrom) {
  pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
  pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
  cumul <- cumul + chrom.lengths[as.character(j)]
  chrom.ends <- c(chrom.ends, cumul)
}
names(chrom.ends) <- names(chrom.lengths)
image(x=c(pos, pos2[length(pos2)]), y=c(pos, pos2[length(pos2)]), z=m, zlim=c(-1,1), col=brewer.pal(10, 'RdBu'), xlab=NA, ylab=NA, main='Correlation Matrix', xaxt='n', yaxt='n')
abline(h=chrom.ends[-length(chrom.ends)], v=chrom.ends[-length(chrom.ends)], lty='dashed')
ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
axis(side=1, at=ax, labels=uni.chrom, lwd=.5, las=1, cex.axis=1)
axis(side=2, at=ax, labels=uni.chrom, lwd=.5, las=1, cex.axis=.5)
axis(side=4, at=ax, labels=uni.chrom, lwd=.5, las=1, cex.axis=.5)
dev.off()

data.info <- dat[, annotationColumns(dat)]

writeData(data.frame(data.info, signif(m, digits=3), check.names=FALSE), file="correlations.tsv")

# EOF
