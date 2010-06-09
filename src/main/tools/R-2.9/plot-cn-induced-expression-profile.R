# ANALYSIS "aCGH tools (beta testing)"/"Plot profiles of matched copy number and expression" (Plot profiles of two priorly matched data sets of copy number and expression. This module must be run on the output from the module Match copy number and expression probes - matched-cn-and-expression.tsv.)
# INPUT GENE_EXPRS matched-cn-and-expression.tsv
# OUTPUT matched-cn-and-expression-profile.png
# PARAMETER sample STRING DEFAULT 1 (The number of the sample to be plotted.)
# PARAMETER chr INTEGER DEFAULT 0 (The chromosome to plot. Use 0 for all.)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# plot-cn-induced-expression-profile.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-05-26

# The plotting command from the intCNGEan package uses par(mfrow) overriding ours.
# Therefore only one profile can be plotted.
# However, the code already contains the necessary functionality to deal with multiple samples,
# in case I find a way how to do it in the future.
# If that happens, the parameter definition will be altered, and this line removed:
samples <- sample

library(CGHcall)
library(intCNGEan)

# read the input file
dat <- read.table('matched-cn-and-expression.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

# check that the input file seems to be coming from the script used to match the two data sets
pos <- c('chromosome','cn.start','cn.end','exp.start','exp.end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This module can only be run on the output file from the module Match copy number and expression probes (matched-cn-and-expression.tsv).')

# build the necessary object
dat$chromosome[dat$chromosome=='X'] <- 23
dat$chromosome[dat$chromosome=='Y'] <- 24
dat$chromosome[dat$chromosome=='MT'] <- 25
dat$chromosome <- as.integer(dat$chromosome)

exprs <- as.matrix(dat[,grep("exprs", names(dat))])
calls <- as.matrix(dat[,grep("flag", names(dat))])
copynumber <- as.matrix(dat[,grep("chip", names(dat))])
segmented <- as.matrix(dat[,grep("segmented", names(dat))])
probloss <- as.matrix(dat[,grep("probloss", names(dat))])
probnorm <- as.matrix(dat[,grep("probnorm", names(dat))])
probgain <- as.matrix(dat[,grep("probgain", names(dat))])

arrays <- sub('exprs\\.', '', colnames(exprs))

cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$cn.start, End=dat$cn.end, row.names=row.names(dat))))
sampleNames(cgh) <- arrays

exp <- new("ExpressionSet", exprs=exprs, featureData=new("AnnotatedDataFrame", data=data.frame(Chromosome=dat$chromosome, Start=dat$exp.start, End=dat$exp.end, row.names=row.names(dat))))
sampleNames(exp) <- arrays

matched <- list(CNdata.matched=cgh, GEdata.matched=exp)

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
to.plot <- to.plot[to.plot<=length(sampleNames(matched$CNdata.matched))]

# check that we have something to plot
if (length(to.plot)==0)
  stop('CHIPSTER-NOTE: Nothing to plot.')

# plot
bitmap(file='matched-cn-and-expression-profile.png', width=image.width/72, height=image.height/72)
if (length(to.plot)==1) {
  intCNGEan.profilesPlot(matched$CNdata.matched, matched$GEdata.matched, sampleNo=to.plot, chr=chr)
} else {
  sq <- sqrt(length(to.plot))
  rows <- ceiling(sq)
  cols <- ceiling(length(to.plot)/rows)
  par(mfrow=c(rows,cols))
  for (sample in to.plot)
    intCNGEan.profilesPlot(matched$CNdata.matched, matched$GEdata.matched, sampleNo=sample, chr=chr)
}
dev.off()

# EOF