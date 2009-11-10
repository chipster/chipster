# ANALYSIS Utilities/"Search by correlation" (Search similarly expressed genes using correlation. 
# Genename should be the row name that appears in all data files, e.g., Affymetrix ID.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT search-correlation.tsv
# PARAMETER genename STRING DEFAULT empty (Gene name)
# PARAMETER correlation.cutoff DECIMAL FROM 0 TO 1 DEFAULT 0.95 (Correlation cut-off for similar expression)
# PARAMETER correlation.method [pearson, spearman] DEFAULT pearson (Correlation method)


# Search genes with correlation to another gene
# JTT 20.6.2006

# Loads the data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Renaming variables
genename<-genename
cutoff<-correlation.cutoff
meth<-correlation.method

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]
gene<-as.numeric(dat2[grep(genename, row.names(dat)),][1,])
len<-dim(dat2)[[1]]
corr<-rep(0, len)
for(i in 1:len) {
   corr[i]<-c(cor(gene, as.numeric(dat2[i,]), method=meth))
}

# Selecting only genes with a correlation coefficient of at least cuttof
dat2<-dat[which(corr<=-cutoff | corr>=cutoff),]
dat2<-data.frame(dat2, cor=corr[which(corr<=-cutoff | corr>=cutoff)])

# Saving the results into a text file
write.table(dat2, "search-correlation.tsv", sep="\t", row.names=T, col.names=T, quote=F)