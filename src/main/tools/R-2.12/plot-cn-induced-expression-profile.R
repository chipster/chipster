# TOOL plot-cn-induced-expression-profile.R: "Plot profiles of matched copy number and expression" (Plot profiles of two priorly matched data sets of copy number and expression. This tool must be run on the output from the tool Match copy number and expression probes - matched-cn-and-expression.tsv.)
# INPUT matched-cn-and-expression.tsv: matched-cn-and-expression.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT matched-cn-and-expression-profile.pdf: matched-cn-and-expression-profile.pdf 
# PARAMETER sample: sample TYPE INTEGER DEFAULT 1 (The number of the sample to be plotted.)
# PARAMETER chromosome: chromosome TYPE INTEGER DEFAULT 0 (The chromosome to plot. Use 0 for all.)
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 6400 DEFAULT 2400 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 6400 DEFAULT 2400 (Height of the plotted network image)

# plot-cn-induced-expression-profile.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2010-10-15

# The plotting command from the intCNGEan package uses par(mfrow) overriding ours.
# Therefore only one profile can be plotted.
# However, the code already contains the necessary functionality to deal with multiple samples,
# in case I find a way how to do it in the future.
# If that happens, the parameter definition will be altered, and this line removed:

library(CGHcall)
library(intCNGEan)

# read the input files
dat <- read.table('matched-cn-and-expression.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
phenodata <- read.table('phenodata.tsv', header=TRUE, sep='\t', as.is=TRUE)

# check if the matched data was produced witha n old version
if (length(grep("^exprs\\.", names(dat)))!=0)
  stop('CHIPSTER-NOTE: The input file matched-cn-and-expression.tsv has been produced with an old version of the Match copy number and expression probes script. Please re-run that script first, and use the output from the new version.')

# check that the input file seems to be coming from the script used to match the two data sets
pos <- c('chromosome','cn.start','cn.end','exp.start','exp.end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This tool can only be run on the output file from the tool Match copy number and expression probes (matched-cn-and-expression.tsv).')

# build the necessary object
dat$chromosome[dat$chromosome=='X'] <- 23
dat$chromosome[dat$chromosome=='Y'] <- 24
dat$chromosome[dat$chromosome=='MT'] <- 25
dat$chromosome <- as.integer(dat$chromosome)

exprs <- as.matrix(dat[,grep("^chip\\.", names(dat))])
calls <- as.matrix(dat[,grep("^flag\\.", names(dat))])
copynumber <- as.matrix(dat[,grep("^copynumber\\.", names(dat))])
segmented <- as.matrix(dat[,grep("^segmented\\.", names(dat))])
probloss <- as.matrix(dat[,grep("^probloss\\.", names(dat))])
probnorm <- as.matrix(dat[,grep("^probnorm\\.", names(dat))])
probgain <- as.matrix(dat[,grep("^probgain\\.", names(dat))])

cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$cn.start, End=dat$cn.end, row.names=row.names(dat))))
sampleNames(cgh) <- phenodata$description_cgh

exp <- new("ExpressionSet", exprs=exprs, featureData=new("AnnotatedDataFrame", data=data.frame(Chromosome=dat$chromosome, Start=dat$exp.start, End=dat$exp.end, row.names=row.names(dat))))
sampleNames(exp) <- phenodata$description

matched <- list(CNdata.matched=cgh, GEdata.matched=exp)

# remove samples that are out of bounds
sample <- sample[sample<=length(sampleNames(matched$CNdata.matched))]

# check that we have something to plot
if (length(sample)==0)
  stop('CHIPSTER-NOTE: Nothing to plot.')

# plot
pdf(file='matched-cn-and-expression-profile.pdf', width=image.width/72, height=image.height/72)
intCNGEan.profilesPlot(matched$CNdata.matched, matched$GEdata.matched, sampleNo=sample, chr=chromosome)
dev.off()

# EOF
