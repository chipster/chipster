# ANALYSIS "aCGH"/"Plot copy number profiles from called aCGH data" (Plot copy number profiles of individual arrays from called aCGH data.)
# INPUT GENE_EXPRS aberrations.tsv, GENERIC phenodata.tsv
# OUTPUT cgh-profile.pdf
# PARAMETER samples STRING DEFAULT 1 (The numbers of the samples to be plotted, separated by commas. Ranges are also supported (e.g. 1,3,7-10).)
# PARAMETER chromosomes STRING DEFAULT 0 (The numbers of the chromosomes to be plotted, separated by commas. 0 means all chromosomes. Ranges are also supported (e.g. 1,3,7-10).)
# PARAMETER image.width INTEGER FROM 200 TO 6400 DEFAULT 2400 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 6400 DEFAULT 2400 (Height of the plotted network image)

# plot-cgh-profile.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-03-28

library(CGHcall)

# read input files
dat <- read.table('aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
phenodata <- read.table("phenodata.tsv", header=TRUE, sep="\t", as.is=TRUE)

# parse samples to be plotted
samples <- gsub('[^0-9,-]', ',', samples)
items <- strsplit(samples, ',')[[1]]
samples.to.plot <- integer()
for (item in items) {
  item <- item[item!='']
  if (length(item)==0) next
  range <- strsplit(item, '-')[[1]]
  range <- range[range!='']
  if (length(range)==0) next
  samples.to.plot <- c(samples.to.plot, seq(range[1], range[length(range)]))
}
samples.to.plot <- unique(samples.to.plot)

# remove samples that are out of bounds
samples.to.plot <- samples.to.plot[samples.to.plot<=nrow(phenodata)]

# check that we have something to plot
if (length(samples.to.plot)==0)
  stop('CHIPSTER-NOTE: Nothing to plot.')

# build a cghCall object from the data
dat$chromosome[dat$chromosome=='X'] <- 23
dat$chromosome[dat$chromosome=='Y'] <- 24
dat$chromosome[dat$chromosome=='MT'] <- 25
dat$chromosome <- as.integer(dat$chromosome)

dat <- dat[dat$chromosome %in% 1:24,]

calls <- as.matrix(dat[,grep("flag", names(dat))])
copynumber <- as.matrix(dat[,grep("chip", names(dat))])
segmented <- as.matrix(dat[,grep("segmented", names(dat))])
probloss <- as.matrix(dat[,grep("probloss", names(dat))])
probnorm <- as.matrix(dat[,grep("probnorm", names(dat))])
probgain <- as.matrix(dat[,grep("probgain", names(dat))])
probamp <- as.matrix(dat[,grep("probamp", names(dat))])

if (ncol(segmented) == 0) {
  cgh <- new('cghRaw', copynumber=copynumber, featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
} else if (ncol(probamp) == 0) {
  cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
} else {
  cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain, probamp=probamp), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
}
sampleNames(cgh) <- phenodata$description

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
  chrs.to.plot <- c(chrs.to.plot, seq(range[1], range[length(range)]))
}
chrs.to.plot <- unique(chrs.to.plot)
chrs.to.plot <- chrs.to.plot[chrs.to.plot %in% dat$chromosome]
if (length(chrs.to.plot)==0)
  chrs.to.plot <- 0

# plot
# pdf(file='cgh-profile.pdf', width=image.width/72, height=image.height/72)
pdf(file='cgh-profile.pdf')
for (sample in samples.to.plot)
  if (0 %in% chrs.to.plot) {
    plot(cgh[,sample]) # dotres=10 -> every 10th log2-ratio is plotted
  } else {
    plot(cgh[chromosomes(cgh) %in% chrs.to.plot, sample]) # dotres=10 -> every 10th log2-ratio is plotted
  }
dev.off()

# EOF
