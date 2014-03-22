# TOOL cna-plot-profile.R: "Plot copy number profiles" (Plot copy number profiles of individual samples.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT cna-profile.pdf: cna-profile.pdf 
# PARAMETER samples: samples TYPE STRING DEFAULT 1 (The numbers of the samples to be plotted, separated by commas. Ranges are also supported (e.g. 1,3,7-10\).)
# PARAMETER chromosomes: chromosomes TYPE STRING DEFAULT 0 (The numbers of the chromosomes to be plotted, separated by commas. 0 means all chromosomes. Ranges are also supported (e.g. 1,3,7-10\).)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-23

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-QDNAseq.R'))

input <- readData("aberrations.tsv")
phenodata <- readPhenodata("phenodata.tsv")

signals <- toQDNAseqSignals(input, chiptype=phenodata$chiptype)
sampleNames(signals) <- phenodata$description

# parse samples to be plotted
if (samples=='0')
  samples <- paste('1-', nrow(phenodata), sep='')
samples <- gsub('[^0-9,-]', ',', samples)
items <- strsplit(samples, ',')[[1]]
samples.to.plot <- integer()
for (item in items) {
  item <- item[item!='']
  if (length(item)==0) next
  range <- strsplit(item, '-')[[1]]
  range <- range[range!='']
  if (length(range)==0) next
  if (length(range) == 1) {
    samples.to.plot <- c(samples.to.plot, as.integer(range))
  } else {
    samples.to.plot <- c(samples.to.plot, seq(as.integer(range[1]), as.integer(range[length(range)])))
  }
}
samples.to.plot <- unique(samples.to.plot)

# remove samples that are out of bounds
samples.to.plot <- samples.to.plot[samples.to.plot <= nrow(phenodata)]

if (0 %in% samples.to.plot)
  samples.to.plot <- 1:nrow(phenodata)

# check that we have something to plot
if (length(samples.to.plot)==0)
  stop('CHIPSTER-NOTE: Nothing to plot.')

# parse chromosomes to be plotted
chromosomes <- gsub('X', '23', chromosomes, ignore.case=TRUE)
chromosomes <- gsub('Y', '24', chromosomes, ignore.case=TRUE)
chromosomes <- gsub('MT', '25', chromosomes, ignore.case=TRUE)
chromosomes <- gsub('[^0-9,-]', ',', chromosomes)
items <- strsplit(chromosomes, ',')[[1]]
chrs.to.plot <- integer()
for (item in items) {
  item <- item[item!='']
  if (length(item)==0) next
  range <- strsplit(item, '-')[[1]]
  range <- range[range!='']
  if (length(range)==0) next
  if (length(range) == 1) {
    chrs.to.plot <- c(chrs.to.plot, as.integer(range))
  } else {
    chrs.to.plot <- c(chrs.to.plot, seq(range[1], range[length(range)]))
  }
}
chrs.to.plot <- unique(chrs.to.plot)
chrs.to.plot <- as.character(chrs.to.plot)
chrs.to.plot[chrs.to.plot == 23] <- 'X'
chrs.to.plot[chrs.to.plot == 24] <- 'Y'
chrs.to.plot[chrs.to.plot == 25] <- 'MT'
chrs.to.plot <- chrs.to.plot[chrs.to.plot %in% input$chromosome]
if (length(chrs.to.plot)==0)
  chrs.to.plot <- 0

# plot
library(png)
tmpfiles <- character(length(samples.to.plot))
for (i in 1:length(samples.to.plot)) {
  tmpfiles[i] <- tempfile()
  bitmap(tmpfiles[i], width=11.7, height=8.3, units='in', res=300)
  if (0 %in% chrs.to.plot) {
    plot(signals[,samples.to.plot[i]])
  } else {
    plot(signals[fData(signals)$chromosome %in% chrs.to.plot, samples.to.plot[i]])
  }
  dev.off()
}

pdf(file='cna-profile.pdf', paper='a4r', width=0, height=0)
par(mai=c(0,0,0,0))
for (i in 1:length(samples.to.plot)) {
  plotPNG = readPNG(tmpfiles[i])
  plot.new()
  rasterImage(plotPNG,0,0,1,1)
}
dev.off()

# EOF
