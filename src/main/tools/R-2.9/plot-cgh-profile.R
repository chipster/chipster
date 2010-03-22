# ANALYSIS "aCGH tools (beta testing)"/"Plot copy number profiles from called aCGH data" (Plot copy number profiles of individual arrays from called aCGH data.)
# INPUT GENE_EXPRS aberrations.tsv, GENERIC phenodata.tsv
# OUTPUT cgh-profile.png
# PARAMETER samples STRING DEFAULT 1 (The numbers of the samples to be plotted, separated by commas. Ranges are also supported (e.g. 1,3,7-10).)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# plot-cgh-profile.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-03-19

library(CGHcall)

# read input files
dat <- read.table('aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")

# parse the input string
samples <- gsub('[^0-9,-]', ',', samples)
items <- strsplit(samples, ',')[[1]]
to.plot <- integer()
for (item in items) {
	item <- item[item!='']
	if (length(item)==0) next
	range <- strsplit(item, '-')[[1]]
	range <- range[range!='']
	if (length(range)==0) next
	to.plot <- c(to.plot, seq(range[1], range[length(range)]))
}
to.plot <- unique(to.plot)

# remove samples that are out of bounds
to.plot <- to.plot[to.plot<=nrow(phenodata)]

# check that we have something to plot
if (length(to.plot)==0)
	stop('Nothing to plot.')

# build a cghCall object from the data
dat$chromosome[dat$chromosome=='X'] <- 23
dat$chromosome[dat$chromosome=='Y'] <- 24
dat$chromosome[dat$chromosome=='MT'] <- 25
dat$chromosome <- as.integer(dat$chromosome)

calls <- as.matrix(dat[,grep("flag", names(dat))])
copynumber <- as.matrix(dat[,grep("chip", names(dat))])
segmented <- as.matrix(dat[,grep("segmented", names(dat))])
probloss <- as.matrix(dat[,grep("probloss", names(dat))])
probnorm <- as.matrix(dat[,grep("probnorm", names(dat))])
probgain <- as.matrix(dat[,grep("probgain", names(dat))])
probamp <- as.matrix(dat[,grep("probamp", names(dat))])

if (ncol(probamp)==0) {
	cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
} else {
	cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain, probamp=probamp), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
}
sampleNames(cgh) <- phenodata$description

# plot
bitmap(file='cgh-profile.png', width=image.width/72, height=image.height/72)
if (length(to.plot)==1) {
	plot(cgh[,to.plot]) # dotres=10 -> every 10th log2-ratio is plotted
} else {
	sq <- sqrt(length(to.plot))
	rows <- ceiling(sq)
	cols <- ceiling(length(to.plot)/rows)
	par(mfrow=c(rows,cols))
	for (sample in to.plot)
		plot(cgh[,sample]) # dotres=10 -> every 10th log2-ratio is plotted
}
dev.off()

# EOF