# TOOL cna-segment-and-call-aberrations.R: "Segment and call copy number aberrations" (This tool segments and calls copy number aberrations from NGS data.)
# INPUT cna-data-table.tsv: "Data table with read counts" TYPE GENE_EXPRS
# OUTPUT aberrations.tsv: aberrations.tsv
# OUTPUT aberration-frequencies.pdf: aberration-frequencies.pdf

# PARAMETER counts: "Counts" TYPE [count: "Original raw counts", corrected: "GC corrected counts"] DEFAULT count (Whether to use original raw read counts, or GC corrected ones.)
# PARAMETER log2: "Log2 transform counts" TYPE [no: no, yes: yes] DEFAULT no (Whether the counts should be log2 transformed.)
# PARAMETER normalization: "Normalization" TYPE [none: none, median: median] DEFAULT none (Normalization method.)
# PARAMETER min.mappability: "Mimimum mappability" TYPE DECIMAL FROM 0 TO 100 DEFAULT 0 (The bins with lower mappability will be removed.)

# add parameters: ploidy, genome.build

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-02-28

source(file.path(chipster.tools.path, 'MPScall', 'CGHcallPlus-R-2.12.R'))
freec.path <- file.path(chipster.tools.path, 'FREEC_Linux64')

input <- read.table('cna-data-table.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
input <- input[order(input$chromosome, input$start),]

ploidy <- 2
binsize <- input$end[1] - input$start[1] + 1
gcfile <- file.path(freec.path, paste('hg19_GC_profile_', binsize / 1000, 'kb.cnp', sep=''))
cat('[general]\n', file='config.txt', sep='')
cat('\n', file='config.txt', sep='', append=TRUE)
cat('ploidy = ', ploidy, '\n', file='config.txt', sep='', append=TRUE)
cat('window = ', binsize, '\n', file='config.txt', sep='', append=TRUE)
cat('breakPointThreshold = 0.6\n', file='config.txt', sep='', append=TRUE)
# cat('telocentromeric = 0\n', file='config.txt', sep='', append=TRUE)
cat('contaminationAdjustment = TRUE\n', file='config.txt', sep='', append=TRUE)
#cat('ontamination = 0.5\n', file='config.txt', sep='', append=TRUE)
# minCNAlength = 1
# minMappabilityPerWindow = 0.85
# uniqueMatch = TRUE
cat('chrLenFile = ', file.path(freec.path, 'hg19.len'), '\n', file='config.txt', sep='', append=TRUE)
cat('gemMappabilityFile = ', file.path(freec.path, 'hg19_GEM_mapp/out50_hg19.gem'), '\n', file='config.txt', sep='', append=TRUE)
if (file.exists(gcfile))
  cat('GCcontentProfile = ', gcfile, '\n', file='config.txt', sep='', append=TRUE)
cat('chrFiles = /home/ischeini/.opt/human/GRCh37\n', file='config.txt', sep='', append=TRUE)
cat('\n', file='config.txt', sep='', append=TRUE)
cat('[sample]\n', file='config.txt', sep='', append=TRUE)
cat('\n', file='config.txt', sep='', append=TRUE)
cat('mateCopyNumberFile = sample.cnp\n', file='config.txt', sep='', append=TRUE)
cat('mateOrientation = 0\n', file='config.txt', sep='', append=TRUE)

anno <- input[,1:3]
ratios <- matrix(nrow=nrow(input), ncol=ncol(input)-3)
colnames(ratios) <- colnames(input)[-(1:3)]
rownames(ratios) <- rownames(input)
calls <- cn <- segmented <- ratios
for (s in colnames(input)[-(1:3)]) {
  cat('Processing: ', s, '\n', sep='')
  write.table(data.frame(input$chromosome, input$start-1, input[,s], stringsAsFactors=FALSE), 'sample.cnp', sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  system(paste(file.path(freec.path, 'freec'), ' -conf config.txt', sep=''))
  out <- read.table('sample.cnp_ratio.txt', header=TRUE, sep='\t', as.is=TRUE)
  ratios[,s] <- out$Ratiofreec <- list(anno=anno, calls=calls, cn=cn, segmented=segmented, ratios=ratios)

  segmented[,s] <- out$MedianRatio
  cn[,s] <- out$CopyNumber
  calls[cn[,s] <  ploidy, s] <- -1
  calls[cn[,s] == ploidy, s] <-  0
  calls[cn[,s] >  ploidy, s] <-  1
  # XX
  # XY
  unlink(c('sample.cnp', 'sample.cnp_CNVs', 'sample.cnp_ratio.txt'))
}
unlink(c('config.txt', 'GC_profile'))
setwd(wd)

anno$chromosome[anno$chromosome == 'X'] <- '23'
anno$chromosome[anno$chromosome == 'Y'] <- '24'
anno$chromosome[anno$chromosome == 'MT'] <- '25'
anno$chromosome <- as.integer(anno$chromosome)

o <- order(anno$chromosome, anno$start)
anno <- anno[o,]
ratios <- as.matrix(ratios[o,])
segmented <- as.matrix(segmented[o,])
cn <- as.matrix(cn[o,])
calls <- as.matrix(calls[o,])
colnames(ratios) <- colnames(segmented) <- colnames(cn) <- colnames(calls) <- colnames(input)[-(1:3)]

ratios[ratios==-1] <- NA
segmented[segmented==-1] <- NA

microarrays <- sprintf('microarray%.3i', 1:ncol(cn))
colnames(ratios) <- paste('chip.', microarrays, sep='')
colnames(segmented) <- paste('segmented.', microarrays, sep='')
colnames(cn) <- paste('copynumber.', microarrays, sep='')
colnames(calls) <- paste('flag.', microarrays, sep='')
dat <- data.frame(bin=paste('bin-', 1:nrow(anno), sep=''), stringsAsFactors=FALSE)
dat <- cbind(dat, anno, calls, cn, ratios, segmented)

frequencyPlot.freec <- 
function (x, main='Frequency Plot', gaincol='blue', losscol='red', misscol=NA, build='GRCh37',... )
{
  chrom <- x$chromosome
  pos <- x$start
  pos2 <- x$end
  chrom.ends <- .getCumulativeChromosomeEnds(build)[1:max(chrom)]
  for (j in 2:max(chrom)) {
    pos[chrom == j] <- pos[chrom == j] + chrom.ends[j-1]
    pos2[chrom == j] <- pos2[chrom == j] + chrom.ends[j-1]
  }
  calls <- as.matrix(dat[,grep("^calls\\.", names(dat))])
  loss.freq <- rowMeans(calls < 0)
  gain.freq <- rowMeans(calls > 0)
  plot(NA, xlim=c(0, max(pos2)), ylim=c(-1,1), type='n', xlab='chromosomes', ylab='frequency', xaxs='i', xaxt='n', yaxs='i', yaxt='n', main=main,...)
  if (!is.na(misscol)) {
    rect(0, -1, max(pos2), 1, col=misscol, border=NA)
    rect(pos, -1, pos2, 1, col='white', border=NA)
  }
  rect(pos, 0, pos2, gain.freq, col=gaincol, border=gaincol)
  rect(pos, 0, pos2, -loss.freq, col=losscol, border=losscol)
  box()
  abline(h=0)
  for (j in 2:max(chrom))
    abline(v=chrom.ends[j-1], lty='dashed')
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1,at=ax,labels=unique(chrom),cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
  mtext('gains', side=2, line=3, at=0.5)
  mtext('losses', side=2, line=3, at=-0.5)
  ### number of data points
  str <- paste(round(nrow(x) / 1000), 'k x ', sep='')
  probe <- median(x$end-x$start+1)
  if (probe < 1000) {
    str <- paste(str, probe, ' bp', sep='')
  } else {
    str <- paste(str, round(probe / 1000), ' kbp', sep='')
  }
  mtext(str, side=3, line=0, adj=0)
}

# write outputs

options(scipen=10)
write.table(dat, file='aberrations.tsv', quote=FALSE, sep='\t', na='', col.names=TRUE, row.names=TRUE)
pdf(file='aberration-frequencies.pdf', paper='a4r', width=0, height=0)
frequencyPlot.freec(dat)
dev.off()

q('no')
# CGHcall

# remove <- rowSums(ratios<0) + rowSums(segmented<0)
remove <- rowSums(is.na(ratios)) + rowSums(is.na(segmented))
anno <- anno[remove==0,]
ratios <- as.matrix(ratios[remove==0,])
segmented <- as.matrix(segmented[remove==0,])
cn <- as.matrix(cn[remove==0,])
calls <- as.matrix(calls[remove==0,])

ratios <- log2(ratios)
segmented <- log2(segmented)

probloss <- calls==-1
probnorm <- calls== 0
probgain <- calls== 1

cgh.raw <- new('cghRaw', assayData=assayDataNew(copynumber=ratios), featureData=new('AnnotatedDataFrame', data=anno))

cgh.seg <- new('cghSeg', assayData=assayDataNew(copynumber=ratios, segmented=segmented), featureData=new('AnnotatedDataFrame', data=anno))

cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=ratios, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=anno))

cgh.all <- cgh
cgh <- cgh[chromosomes(cgh) < 23,]

# EOF
