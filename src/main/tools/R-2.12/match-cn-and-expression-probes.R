# TOOL match-cn-and-expression-probes.R: "Match copy number and expression probes" (Matches the probes of a copy number data set with probes of an expression data set, using their chromosomal locations. Running this tool is a prerequisite for testing copy-number-induced effects on expression.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT phenodata_cgh.tsv: phenodata_cgh.tsv TYPE GENERIC 
# INPUT phenodata_exp.tsv: phenodata_exp.tsv TYPE GENERIC 
# OUTPUT matched-cn-and-expression.tsv: matched-cn-and-expression.tsv 
# OUTPUT matched-cn-and-expression-heatmap.pdf: matched-cn-and-expression-heatmap.pdf 
# OUTPUT matched-phenodata.tsv: matched-phenodata.tsv 
# PARAMETER sample.identifiers.1: sample.identifiers.1 TYPE METACOLUMN_SEL DEFAULT Sample (The phenodata column for data set 1 used to link the two data sets together.)
# PARAMETER sample.identifiers.2: sample.identifiers.2 TYPE METACOLUMN_SEL DEFAULT Sample (The phenodata column for data set 2 used to link the two data sets together.)
# PARAMETER method: method TYPE [distance: distance, overlap: overlap, overlapplus: overlapplus] DEFAULT distance (The method for linking copy number and expression probes together.)

# match-cn-and-expression-probes.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-03-28

library(CGHcall)
library(intCNGEan)

dat1 <- read.table('aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
dat2 <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
phenodata1 <- read.table("phenodata_cgh.tsv", header=T, sep='\t', as.is=TRUE)
phenodata2 <- read.table("phenodata_exp.tsv", header=T, sep='\t', as.is=TRUE)

# determine which dataset is cgh
if (length(grep("^probnorm\\.", names(dat1)))!=0) {
  cgh <- dat1
  exp <- dat2
  phenodata_cgh <- phenodata1
  phenodata_exp <- phenodata2
  samples_cgh <- sample.identifiers.1
  samples_exp <- sample.identifiers.2
} else if (length(grep("^probnorm\\.", names(dat2)))!=0) {
  cgh <- dat2
  exp <- dat1
  phenodata_cgh <- phenodata2
  phenodata_exp <- phenodata1
  samples_cgh <- sample.identifiers.2
  samples_exp <- sample.identifiers.1
} else {
  stop('CHIPSTER-NOTE: Could not detect the aCGH data set.')
}

# check that both data sets have the probe position information
pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(cgh)))!=0)
  stop('CHIPSTER-NOTE: Both data sets must have the following columns: chromosome, start, end.')
if (length(setdiff(pos, colnames(exp)))!=0)
  stop('CHIPSTER-NOTE: Both data sets must have the following columns: chromosome, start, end.')

# check for unambiguity of sample identifiers
if (nrow(phenodata_cgh)!=length(unique(phenodata_cgh[,samples_cgh])))
  stop('CHIPSTER-NOTE: Unambigous aCGH sample identifiers: ', paste(phenodata_cgh[duplicated(phenodata_cgh[,samples_cgh]),samples_cgh], collapse=', ')) 
if (nrow(phenodata_exp)!=length(unique(phenodata_exp[,samples_exp])))
  stop('CHIPSTER-NOTE: Unambigous expression sample identifiers: ', paste(phenodata_exp[duplicated(phenodata_exp[,samples_exp]),samples_exp], collapse=', ')) 

common.samples <- intersect(phenodata_cgh[,samples_cgh], phenodata_exp[,samples_exp])
rownames(phenodata_cgh) <- phenodata_cgh[,samples_cgh]
rownames(phenodata_exp) <- phenodata_exp[,samples_exp]
phenodata_cgh$n <- 1:nrow(phenodata_cgh)
phenodata_exp$n <- 1:nrow(phenodata_exp)
index_cgh <- phenodata_cgh[common.samples, 'n']
index_exp <- phenodata_exp[common.samples, 'n']

# build a cghCall object from the cgh data
cgh$chromosome[cgh$chromosome=='X'] <- 23
cgh$chromosome[cgh$chromosome=='Y'] <- 24
cgh$chromosome[cgh$chromosome=='MT'] <- 25
cgh$chromosome <- as.integer(cgh$chromosome)
cgh <- cgh[!is.na(cgh$chromosome),]
calls <- as.matrix(cgh[,grep("^flag\\.", names(cgh))])[,index_cgh]
copynumber <- as.matrix(cgh[,grep("^chip\\.", names(cgh))])[,index_cgh]
segmented <- as.matrix(cgh[,grep("^segmented\\.", names(cgh))])[,index_cgh]
probloss <- as.matrix(cgh[,grep("^probloss\\.", names(cgh))])[,index_cgh]
probnorm <- as.matrix(cgh[,grep("^probnorm\\.", names(cgh))])[,index_cgh]
probgain <- as.matrix(cgh[,grep("^probgain\\.", names(cgh))])[,index_cgh]
if (2 %in% calls) {
  probamp <- as.matrix(cgh[,grep("^probamp\\.", names(cgh))])[,index_cgh]
  probgain <- probgain + probamp
  calls[calls==2] <- 1
}
cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=cgh$chromosome, Start=cgh$start, End=cgh$end, row.names=row.names(cgh))))
sampleNames(cgh) <- phenodata_cgh[common.samples, samples_cgh]

# build an ExpressionSet object from the expression data
exp$chromosome[exp$chromosome=='X'] <- 23
exp$chromosome[exp$chromosome=='Y'] <- 24
exp$chromosome[exp$chromosome=='MT'] <- 25
exp$chromosome <- as.integer(exp$chromosome)
exp <- exp[!is.na(exp$chromosome),]
exprs <- as.matrix(exp[,grep("^chip\\.", names(exp))])[,index_exp]
exp <- new("ExpressionSet", exprs=exprs, featureData=new("AnnotatedDataFrame", data=data.frame(Chromosome=exp$chromosome, Start=exp$start, End=exp$end, row.names=row.names(exp))))
sampleNames(exp) <- phenodata_exp[common.samples, samples_exp]

# match probes
matched <- intCNGEan.match(cgh, exp, CNbpend='yes', GEbpend='yes', method=method)

# plot heatmaps
pdf(file='matched-cn-and-expression-heatmap.pdf')
intCNGEan.heatmaps(matched$CNdata.matched, matched$GEdata.matched)
dev.off()

# separate elements from the resulting object

dat3 <- data.frame(matched$CNdata.matched@featureData@data, matched$GEdata.matched@featureData@data)
colnames(dat3) <- c('chromosome', 'cn.start', 'cn.end', 'exp.probe','exp.start', 'exp.end')
dat3$exp.probe <- rownames(matched$GEdata.matched@featureData@data)

dat3$loss.freq <- round(mean(as.data.frame(t(assayDataElement(matched$CNdata.matched, "calls")==-1))), digits=3)
dat3$gain.freq <- round(mean(as.data.frame(t(assayDataElement(matched$CNdata.matched, "calls")==1))), digits=3)

exprs <- assayDataElement(matched$GEdata.matched, 'exprs')
samples <- microarrays <- sprintf('microarray%.3i', 1:ncol(exprs))
colnames(exprs) <- paste('chip.', samples, sep='')
dat3 <- cbind(dat3, exprs)

calls <- assayDataElement(matched$CNdata.matched, 'calls')
colnames(calls) <- paste('flag.', samples, sep='')
dat3 <- cbind(dat3, calls)

copynumber <- assayDataElement(matched$CNdata.matched, 'copynumber')
colnames(copynumber) <- paste('copynumber.', samples, sep='')
dat3 <- cbind(dat3, copynumber)

segmented <- assayDataElement(matched$CNdata.matched, 'segmented')
colnames(segmented) <- paste('segmented.', samples, sep='')
dat3 <- cbind(dat3, segmented)

probloss <- assayDataElement(matched$CNdata.matched, 'probloss')
colnames(probloss) <- paste('probloss.', samples, sep='')
dat3 <- cbind(dat3, probloss)

probnorm <- assayDataElement(matched$CNdata.matched, 'probnorm')
colnames(probnorm) <- paste('probnorm.', samples, sep='')
dat3 <- cbind(dat3, probnorm)

probgain <- assayDataElement(matched$CNdata.matched, 'probgain')
colnames(probgain) <- paste('probgain.', samples, sep='')
dat3 <- cbind(dat3, probgain)

dat3$chromosome <- as.character(dat3$chromosome)
dat3$chromosome[dat3$chromosome=='23'] <- 'X'
dat3$chromosome[dat3$chromosome=='24'] <- 'Y'
dat3$chromosome[dat3$chromosome=='25'] <- 'MT'

# add additional phenodata columns from phenodata_cgh to phenodata_exp
phenodata_cgh <- phenodata_cgh[index_cgh,]
phenodata_exp <- phenodata_exp[index_exp,]
phenodata_exp$description_cgh <- phenodata_cgh$description
for (col in setdiff(colnames(phenodata_cgh), colnames(phenodata_exp)))
  phenodata_exp[,col] <- phenodata_cgh[,col]
phenodata_exp$sample <- samples
phenodata_exp$n <- NULL

# write output
write.table(dat3, file='matched-cn-and-expression.tsv', quote=FALSE, sep='\t')
write.table(phenodata_exp, file='matched-phenodata.tsv', quote=FALSE, sep='\t', na='', row.names=FALSE)

# EOF
