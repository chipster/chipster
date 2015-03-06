# TOOL acgh-plot-profile.R: "Plot copy number profiles" (Plot copy number profiles of individual samples.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT cgh-profile.pdf: cgh-profile.pdf 
# PARAMETER samples: Samples TYPE STRING DEFAULT 1 (The numbers of the samples to be plotted, separated by commas. Ranges are also supported (e.g. 1,3,7-10\).)
# PARAMETER chromosomes: Chromosomes TYPE STRING DEFAULT 0 (The numbers of the chromosomes to be plotted, separated by commas. 0 means all chromosomes. Ranges are also supported (e.g. 1,3,7-10\).)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-23

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-CGHcall.R'))

input <- readData("aberrations.tsv")
phenodata <- readPhenodata("phenodata.tsv")

cgh <- toCgh(input)
sampleNames(cgh) <- phenodata$description

samples.to.plot <- parseSamplesToPlot(samples, 1:nrow(phenodata))
chrs.to.plot <- parseChromosomesToPlot(chromosomes, chromosomes(cgh))

# plot
library(png)
tmpfiles <- character(length(samples.to.plot))
for (i in 1:length(samples.to.plot)) {
  tmpfiles[i] <- tempfile()
  bitmap(tmpfiles[i], width=11.7, height=8.3, units='in', res=300)
  plot(cgh[chrs.to.plot, samples.to.plot[i]])
  dev.off()
}

pdf(file='cgh-profile.pdf', paper='a4r', width=0, height=0)
par(mai=c(0,0,0,0))
for (i in 1:length(samples.to.plot)) {
  plotPNG = readPNG(tmpfiles[i])
  plot.new()
  rasterImage(plotPNG,0,0,1,1)
}
dev.off()

# EOF
