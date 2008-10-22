# ANALYSIS Statistics/"Correlate with phenodata" (Finds genes that correlate highly with the specified phenodata column.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT correlate-with-phenodata.tsv
# PARAMETER correlation.cutoff DECIMAL FROM 0 TO 1 DEFAULT 0.95 (Correlation cut-off for similar expression)
# PARAMETER correlation.method [pearson, spearman] DEFAULT pearson (Correlation method)
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column containing the data to correlate genes with)


# Find genes that correlate with phenodata
# JTT 19.7.2007

# Renaming variables
cutoff<-correlation.cutoff
meth<-correlation.method

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Test needs a parameter "groups" that specifies the grouping of the samples
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
groups<-phenodata[,grep(column, colnames(phenodata))]

# Sanity checks
if(nrow(phenodata)!=ncol(dat2)) {
   stop("Dimensions of data and phenodata do not match! Can't execute the analysis.")
}

# Calculating the correlations
p<-rep(0, nrow(dat2))
for(i in 1:nrow(dat2)) {
   p[i]<-cor(groups, as.numeric(dat2[i,]), method=meth)
}

# Filtering the data on the correlation
dat<-dat[abs(p)>=cutoff,]
p<-p[abs(p)>=cutoff]

# Saving the results
write.table(data.frame(dat, corr=round(p, digits=4)), file="correlate-with-phenodata.tsv", sep="\t", row.names=T, col.names=T, quote=F)


