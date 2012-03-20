# TOOL cna-segment-and-call-aberrations.R: "Segment and call copy number aberrations" (This tool segments and calls copy number aberrations from NGS data.)
# INPUT cna-data-table.tsv: "Data table with read counts" TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT aberrations.tsv: aberrations.tsv
# OUTPUT aberration-frequencies.pdf: aberration-frequencies.pdf
# PARAMETER organism: "Organism" TYPE [human: human] DEFAULT human (Organism.)
# PARAMETER genome.build: "Genome build" TYPE [GRCh37: GRCh37] DEFAULT GRCh37 (Genome build.)
# PARAMETER breakpoint.threshold: "Breakpoint threshold" TYPE DECIMAL DEFAULT 0.8 (Threshold for segmentation of normalized profiles. Use something like 0.6 to get more segments, and thus more predicted CNVs.)
# PARAMETER min.mappability: "Mimimum mappability" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.85 (The bins with lower mappability will be removed.)
# PARAMETER ploidy: "Phenodata column with ploidy of samples" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column with ploidy)
# PARAMETER cellularity: "Phenodata column with cellularity of samples" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column with the proportion of tumor cells in the sample)
# PARAMETER unique.match: "Unique match" TYPE [FALSE: no, TRUE: yes] DEFAULT FALSE (Use a mappability profile to correct read counts.)

# add parameters: ploidy, read length, breakPointThreshold

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-03-20

source(file.path(chipster.tools.path, 'MPScall', 'CGHcallPlus-R-2.12.R'))
freec.path <- file.path(chipster.tools.path, 'FREEC_Linux64')
options(scipen=10)

input <- read.table('cna-data-table.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
input <- input[input$chromosome %in% c(1:22, 'X', 'Y'),]
input <- input[order(input$chromosome, input$start),]
binsize <- input$end[1] - input$start[1] + 1

phenodata <- read.table('phenodata.tsv', header=TRUE, sep='\t')

if (ploidy == 'EMPTY') {
  ploidy <- rep(2, ncol(input)-3)
} else {
  ploidy <- as.numeric(phenodata[,ploidy])
  ploidy[is.na(ploidy)] <- 2
}
names(ploidy) <- colnames(input)[-(1:3)]
if (cellularity == 'EMPTY') {
  cellularity <- rep(1, ncol(input)-3)
} else {
  cellularity <- as.numeric(phenodata[,cellularity])
  cellularity[is.na(cellularity)] <- 1
}
names(cellularity) <- colnames(input)[-(1:3)]
contamination <- 1 - cellularity

anno <- input[,1:3]
ratios <- matrix(nrow=nrow(input), ncol=ncol(input)-3)
colnames(ratios) <- colnames(input)[-(1:3)]
rownames(ratios) <- rownames(input)
calls <- cn <- segmented <- ratios
for (s in colnames(input)[-(1:3)]) {
  cat('[general]\n', file='config.txt', sep='')
  cat('\n', file='config.txt', sep='', append=TRUE)
  cat('breakPointThreshold = ', breakpoint.threshold, '\n', file='config.txt', sep='', append=TRUE)
  cat('chrLenFile = ', file.path(freec.path, genome.build, 'chr.len'), '\n', file='config.txt', sep='', append=TRUE)
  if (contamination[s] > 0) {
    cat('contamination = ', contamination[s], '\n', file='config.txt', sep='', append=TRUE)
    cat('contaminationAdjustment = TRUE\n', file='config.txt', sep='', append=TRUE)
  }
  cat('GCcontentProfile = ', file.path(freec.path, genome.build, paste('GC_profile_', binsize / 1000, 'kb.cnp', sep='')), '\n', file='config.txt', sep='', append=TRUE)
  cat('gemMappabilityFile = ', file.path(freec.path, genome.build, 'GEM_mapp', 'out50.gem'), '\n', file='config.txt', sep='', append=TRUE)
  # cat('minCNAlength = 1\n', file='config.txt', sep='', append=TRUE)
  cat('minMappabilityPerWindow = ', min.mappability, '\n', file='config.txt', sep='', append=TRUE)
  cat('ploidy = ', ploidy[s], '\n', file='config.txt', sep='', append=TRUE)
  # cat('telocentromeric = 50000\n', file='config.txt', sep='', append=TRUE)
  cat('uniqueMatch = ', unique.match, '\n', file='config.txt', sep='', append=TRUE)
  cat('window = ', binsize, '\n', file='config.txt', sep='', append=TRUE)
  cat('\n', file='config.txt', sep='', append=TRUE)
  cat('[sample]\n', file='config.txt', sep='', append=TRUE)
  cat('\n', file='config.txt', sep='', append=TRUE)
  cat('mateCopyNumberFile = sample.cnp\n', file='config.txt', sep='', append=TRUE)
  cat('mateOrientation = 0\n', file='config.txt', sep='', append=TRUE)

  write.table(data.frame(input$chromosome, input$start-1, input[,s], stringsAsFactors=FALSE), 'sample.cnp', sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  system(paste(file.path(freec.path, 'freec'), ' -conf config.txt', sep=''))
  out <- read.table('sample.cnp_ratio.txt', header=TRUE, sep='\t', as.is=TRUE)
  ratios[,s] <- out$Ratio
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

ratios[ratios == -1] <- NA
segmented[segmented == -1] <- NA
ratios <- log2(ratios)
segmented <- log2(segmented)
ratios[ratios == -10] <- NA # or NA
segmented[segmented == -10] <- NA # or NA

identifiers <- sub('^chip\\.', '', colnames(ratios))
colnames(segmented) <- paste('segmented.', identifiers, sep='')
colnames(cn) <- paste('copynumber.', identifiers, sep='')
colnames(calls) <- paste('flag.', identifiers, sep='')
dat <- cbind(anno, calls, cn, ratios, segmented)

.getChromosomeLengths <- function(build) {
    build <- as.integer(gsub('[^0-9]', '', build))
    if (build == 34 || build == 16) {
       chromosome.lengths <- c(246127941, 243615958, 199344050, 191731959, 181034922, 170914576, 158545518, 146308819, 136372045, 135037215, 134482954, 132078379, 113042980, 105311216, 100256656, 90041932, 81860266, 76115139, 63811651, 63741868, 46976097, 49396972, 153692391, 50286555)
    } else if (build == 35 || build == 17) {
       chromosome.lengths <- c(245522847, 243018229, 199505740, 191411218, 180857866, 170975699, 158628139, 146274826, 138429268, 135413628, 134452384, 132449811, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49554710, 154824264, 57701691)
    } else if (build == 36 || build == 18) {
       chromosome.lengths <- c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
    } else {
       chromosome.lengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
    }
    names(chromosome.lengths) <- 1:24
    chromosome.lengths
}

frequencyPlot.freec <- 
function (x, main='Frequency Plot', gaincol='blue', losscol='red', misscol=NA, build='GRCh37',... )
{
  chrom <- x$chromosome
  pos <- x$start
  pos2 <- x$end
  uni.chrom <- unique(chrom)
  chrom.lengths <- .getChromosomeLengths(build)[as.character(uni.chrom)]
  chrom.ends <- integer()
  cumul <- 0
  for (j in uni.chrom) {
    pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
    pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
    cumul <- cumul + chrom.lengths[as.character(j)]
    chrom.ends <- c(chrom.ends, cumul)
  }
  names(chrom.ends) <- names(chrom.lengths)
  calls <- as.matrix(dat[,grep("^flag\\.", names(dat))])
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
  if (length(chrom.ends) > 1)
    for (j in names(chrom.ends)[-length(chrom.ends)])
      abline(v=chrom.ends[j], lty='dashed')
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1,at=ax,labels=uni.chrom,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
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

pdf(file='aberration-frequencies.pdf', paper='a4r', width=0, height=0)
frequencyPlot.freec(dat)
dev.off()

dat$chromosome <- as.character(dat$chromosome)
dat$chromosome[dat$chromosome=='23'] <- 'X'
dat$chromosome[dat$chromosome=='24'] <- 'Y'
dat$chromosome[dat$chromosome=='25'] <- 'MT'

write.table(dat, file='aberrations.tsv', quote=FALSE, sep='\t')

# EOF
