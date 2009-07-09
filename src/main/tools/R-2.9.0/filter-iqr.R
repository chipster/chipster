# ANALYSIS Preprocessing/"Filter by interquartile range" (Filter genes by their interquartile range.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT iqr-filter.tsv
# PARAMETER iqr.threshold DECIMAL FROM 0 TO 10 DEFAULT 0.5 (Interquartile range)


# JTT, 24.10.2007

# Renaming variables
iqrt<-iqr.threshold

# Loading the libraries
library(genefilter)

# Loads the normalized data
# This needs a parameter "file" that tells which data file should be loaded.
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]
calls<-data.frame(calls)
dat2<-data.frame(dat2)

# Sanity checks
if(ncol(dat2)==1) {
   stop("You need to have at least two chips to apply filtering by genes!")
}

# Filter

f<-function(x) {
   (IQR(x)>=iqrt)
}
ff<-filterfun(f)
sel<-genefilter(dat2, ff)
set<-dat[sel, ]

# Saving the results
write.table(data.frame(set), file=("iqr-filter.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
