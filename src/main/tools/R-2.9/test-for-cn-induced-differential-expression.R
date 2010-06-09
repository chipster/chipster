# ANALYSIS "aCGH tools (beta testing)"/"Test for copy number induced expression changes" (Nonparametric testing for changes in expression induced by a change in DNA copy number. The copy number and expression probes of the two data sets must be matched together beforehand using the Match copy number and expression probes module.)
# INPUT GENE_EXPRS matched-cn-and-expression.tsv
# OUTPUT cn-induced-expression.tsv
# PARAMETER test.statistic [wcvm, wmw] DEFAULT wcvm (The test statistic to use.)
# PARAMETER analysis.type [univariate, regional] DEFAULT univariate (The type of the analysis.)
# PARAMETER nperm INTEGER DEFAULT 10000 (The number of permutations used for the p-value calculation.)

# test-for-cn-induced-differential-expression.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-06-04

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
rownames(exprs) <- dat$exp.probe
arrays <- sub('exprs\\.', '', colnames(exprs))

calls <- as.matrix(dat[,grep("flag", names(dat))])
copynumber <- as.matrix(dat[,grep("chip", names(dat))])
segmented <- as.matrix(dat[,grep("segmented", names(dat))])
probloss <- as.matrix(dat[,grep("probloss", names(dat))])
probnorm <- as.matrix(dat[,grep("probnorm", names(dat))])
probgain <- as.matrix(dat[,grep("probgain", names(dat))])

cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$cn.start, End=dat$cn.end, row.names=row.names(dat))))
sampleNames(cgh) <- arrays

exp <- new("ExpressionSet", exprs=exprs, featureData=new("AnnotatedDataFrame", data=data.frame(Chromosome=dat$chromosome, Start=dat$exp.start, End=dat$exp.end, row.names=dat$exp.probe)))
sampleNames(exp) <- arrays

matched <- list(CNdata.matched=cgh, GEdata.matched=exp)

# tune and test
tuned <- intCNGEan.tune(matched$CNdata.matched, matched$GEdata.matched, test.statistic=test.statistic)
result <- intCNGEan.test(tuned, analysis.type=analysis.type, test.statistic=test.statistic, nperm=nperm)

if (any(rownames(tuned$ann) != rownames(result)))
  stop('CHIPSTER-NOTE: The rownames for the objects tuned and result do not match. Please report this to Ilari Scheinin.')

# format and write result table
colnames(result)[1:3] <- tolower(colnames(result)[1:3])
result$probes <- rownames(tuned$datafortest)
arrays <- colnames(tuned$datafortest)[(2*tuned$nosamp+1):(3*tuned$nosamp)]
colnames(tuned$datafortest) <- c(paste('prob.1.', arrays, sep=''), paste('prob.2.', arrays, sep=''), paste('chip.', arrays, sep=''))
result <- cbind(result, tuned$datafortest)
result <- result[order(result$adj.p),]

result$chromosome <- as.character(result$chromosome)
result$chromosome[result$chromosome=='23'] <- 'X'
result$chromosome[result$chromosome=='24'] <- 'Y'
result$chromosome[result$chromosome=='25'] <- 'MT'

write.table(result, file='cn-induced-expression.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF