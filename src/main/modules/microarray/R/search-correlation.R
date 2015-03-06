# TOOL search-correlation.R: "Search by correlation" (Search similarly expressed genes using correlation. Genename should be the row name that appears in all data files, e.g., Affymetrix ID.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT search-correlation.tsv: search-correlation.tsv 
# PARAMETER genename: "Gene name" TYPE STRING DEFAULT empty (Gene name)
# PARAMETER correlation.cutoff: "Correlation cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.95 (Correlation cut-off for similar expression)
# PARAMETER correlation.method: "Correlation method" TYPE [pearson: Pearson, spearman: Spearman] DEFAULT pearson (Correlation method)


# JTT 20.6.2006: Search genes with correlation to another gene
# MK 25.10.2013: Modified to search also symbols

# Parameter settings (default) for testing purposes
#correlation.cutoff<-0.95
#correlation.method<-"pearson"
#genename<-"1007_s_at"

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

gene <- NULL;
if(length(grep(genename, row.names(dat)))>=1) {
	gene<-as.numeric(dat2[grep(genename, row.names(dat)),][1,])	
} else {
	if(length(grep("symbol", colnames(dat)))==1) {
		if(length(grep(genename, dat[,"symbol"]))>=1) {
			gene<-as.numeric(dat2[grep(genename, dat[,"symbol"]),][1,])
		}
	}
}

if(is.null(gene)) {
	stop("CHIPSTER-NOTE: Search string not found or found too many times")
}

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
