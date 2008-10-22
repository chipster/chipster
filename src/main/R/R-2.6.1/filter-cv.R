# ANALYSIS Preprocessing/"Filter by CV" (Filter genes by their coefficient of variation. The filtering is automatically
# performed using the median of the CV values.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT cv-filter.tsv


# JTT, 14.1.20078

# Loading the libraries
library(genefilter)

# Loads the normalized data
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
colno<-ncol(dat2)
g.mean <- rowSums(dat2)/colno
g.sd <- rowSds(dat2)
g.cv <- g.sd / g.mean
ffun <- filterfun( cv(0,median(g.cv) ))
sel  <- genefilter(dat2, ffun)
set <- dat[sel, ]

# Saving the results
write.table(data.frame(set), file=("cv-filter.tsv"), sep="\t", row.names=T, col.names=T, quote=F)

