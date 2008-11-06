# ANALYSIS Utilities/"Calculate descriptive statistics" (Calculates basic descriptive statistics for all genes.
# These include parametric and non-parametric location and spread descriptives.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT descr-stats.tsv
# PARAMETER calculate.descriptives.for [genes, chips] DEFAULT genes (Descriptive statistics are calculated for...)


# Two-group parametric and non-parametric tests
# JTT 31.1.2008

# Loads the libraries
library(genefilter)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Descriptives for chips
if(calculate.descriptives.for=="chips") {
   dat2<-t(dat2)
}

# Calculates the descriptive statistics
nc<-ncol(dat2)
rsum<-rowSums(dat2)
rmean<-rsum/nc
rsd<-rowSds(dat2)
rcv<-rsd/rmean
rstexp<-rmean/rsd
rmin<-apply(X=dat2, MARGIN=1, FUN=min)
rmax<-apply(X=dat2, MARGIN=1, FUN=max)
rrange<-rmax-rmin
rmedian<-apply(X=dat2, MARGIN=1, FUN=median)
riqr<-apply(X=dat2, MARGIN=1, FUN=IQR)

# Saving the results
if(calculate.descriptives.for=="chips") {
   write.table(data.frame(rownames(dat2), average=rmean, median=rmedian, sd=rsd, cv=rcv, st.exp=rstexp, min=rmin, max=rmax, range=rrange, iqr=riqr), file="descr-stats.tsv", sep="\t", row.names=T, col.names=T, quote=F)
} else {
   write.table(data.frame(dat, average=rmean, median=rmedian, sd=rsd, cv=rcv, st.exp=rstexp, min=rmin, max=rmax, range=rrange, iqr=riqr), file="descr-stats.tsv", sep="\t", row.names=T, col.names=T, quote=F)
   write.table(data.frame(chip.average=rmean, chip.median=rmedian, chip.sd=rsd, chip.cv=rcv, chip.st.exp=rstexp, chip.min=rmin, chip.max=rmax, chip.range=rrange, chip.iqr=riqr), file="descriptives.tsv", sep="\t", row.names=T, col.names=T, quote=F)

}
