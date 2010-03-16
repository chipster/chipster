# ANALYSIS "aCGH tools (beta testing)"/"Test for DNA copy number induced differential expression" (Test for changes in expression induced by a change in DNA copy number.)
# INPUT GENE_EXPRS aberrations.tsv, GENE_EXPRS normalized.tsv, GENERIC phenodata_cgh.tsv, GENERIC phenodata_exp.tsv
# OUTPUT cn-induced-expression.tsv, cn-induced-expression-heatmap.png, cn-induced-expression.RData
# PARAMETER samples1 METACOLUMN_SEL DEFAULT EMPTY (The phenodata column for data set 1 used to link the two data sets together.)
# PARAMETER samples2 METACOLUMN_SEL DEFAULT EMPTY (The phenodata column for data set 2 used to link the two data sets together.)
# PARAMETER method [distance, overlap, overlapplus] DEFAULT distance (The method for linking copy number and expression probes together.)
# PARAMETER test.statistic [wcvm, wmw] DEFAULT wcvm (The test statistic to use.)
# PARAMETER analysis.type [univariate, regional] DEFAULT univariate (The type of the analysis.)

# test-for-cn-induced-differential-expression.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-03-11

library(CGHcall)
library(intCNGEan)

# dat1 <- read.table('aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
# dat2 <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
dat1 <- read.table('aberrations.tsv', header=TRUE, sep='\t', row.names=1)
dat2 <- read.table('normalized.tsv', header=TRUE, sep='\t', row.names=1)
phenodata1 <- read.table("phenodata_cgh.tsv", header=T, sep="\t")
phenodata2 <- read.table("phenodata_exp.tsv", header=T, sep="\t")

# determine which dataset is cgh
if (length(grep("probnorm", names(dat1)))!=0) {
  cgh <- dat1
  exp <- dat2
  phenodata_cgh <- phenodata1
  phenodata_exp <- phenodata2
  samples_cgh <- samples1
  samples_exp <- samples2
} else if (length(grep("probnorm", names(dat2)))!=0) {
  cgh <- dat2
  exp <- dat1
  phenodata_cgh <- phenodata2
  phenodata_exp <- phenodata1
  samples_cgh <- samples2
  samples_exp <- samples1
} else {
  stop('Could not detect the aCGH data set.')
}

# check for unambiguity of sample identifiers
if (nrow(phenodata_cgh)!=length(unique(phenodata_cgh[,samples_cgh])))
  stop('Unambigous sample identifiers: ', paste(phenodata_cgh[,samples_cgh], collapse=', ')) 
if (nrow(phenodata_exp)!=length(unique(phenodata_exp[,samples_exp])))
  stop('Unambigous sample identifiers: ', paste(phenodata_exp[,samples_exp], collapse=', ')) 

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
calls <- as.matrix(cgh[,grep("flag", names(cgh))])[,index_cgh]
copynumber <- as.matrix(cgh[,grep("chip", names(cgh))])[,index_cgh]
segmented <- as.matrix(cgh[,grep("segmented", names(cgh))])[,index_cgh]
probloss <- as.matrix(cgh[,grep("probloss", names(cgh))])[,index_cgh]
probnorm <- as.matrix(cgh[,grep("probnorm", names(cgh))])[,index_cgh]
probgain <- as.matrix(cgh[,grep("probgain", names(cgh))])[,index_cgh]
probamp <- as.matrix(cgh[,grep("probamp", names(cgh))])
if (ncol(probamp)==0) {
  cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=cgh$chromosome, Start=cgh$start, End=cgh$end, row.names=row.names(cgh))))
} else {
  probamp <- probamp[,index_cgh]
  cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain, probamp=probamp), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=cgh$chromosome, Start=cgh$start, End=cgh$end, row.names=row.names(cgh))))
}
sampleNames(cgh) <- phenodata_cgh[common.samples, samples_cgh]

# build an ExpressionSet object from the expression data
exp$chromosome[exp$chromosome=='X'] <- 23
exp$chromosome[exp$chromosome=='Y'] <- 24
exp$chromosome[exp$chromosome=='MT'] <- 25
exp$chromosome <- as.integer(exp$chromosome)
exp <- exp[!is.na(exp$chromosome),]
exprs <- as.matrix(exp[,grep("chip", names(exp))])[,index_exp]
exp <- new("ExpressionSet", exprs=exprs, featureData=new("AnnotatedDataFrame", data=data.frame(Chromosome=exp$chromosome, Start=exp$start, End=exp$end, row.names=row.names(exp))))
sampleNames(exp) <- phenodata_exp[common.samples, samples_exp]

# test for copy-number-induced effects on expression
matched <- intCNGEan.match(cgh, exp, CNbpend='yes', GEbpend='yes', method='distance')
tuned <- intCNGEan.tune(matched$CNdata.matched, matched$GEdata.matched, test.statistic=test.statistic)
result <- intCNGEan.test(tuned, analysis.type=analysis.type, test.statistic=test.statistic)
colnames(result)[1:3] <- tolower(colnames(result)[1:3])

# ordering the results might mess up the plotting functions, so it has been commented out
# result <- result[order(result$adj.p),]

write.table(result, file='cn-induced-expression.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# there seems to be a bug in the package, which causes the plot to be two pages, with the real stuff on the second one.
# so we'll plot it in two files
bitmap(file='cn-induced-expression-heatmap%d.png', width=600/72, height=600/72)
intCNGEan.heatmaps(matched$CNdata.matched, matched$GEdata.matched)
dev.off()

# in case the bug is fixed in the future, let's first check if we indeed have two files.
# if yes, let's take the second one.
# if not, let's just take the first and only file.
if (file.exists('cn-induced-expression-heatmap2.png')) {
  file.rename('cn-induced-expression-heatmap2.png', 'cn-induced-expression-heatmap.png')
  file.remove('cn-induced-expression-heatmap1.png')
} else {
  file.rename('cn-induced-expression-heatmap1.png', 'cn-induced-expression-heatmap.png')
}

save(matched, tuned, result, file='cn-induced-expression.RData')

# EOF