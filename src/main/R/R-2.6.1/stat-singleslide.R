# ANALYSIS Statistics/"Single-slide methods" (These are applicable if you only have a single cDNA slide. The noise-envelope
# uses SD for filtering the genes. Newton's method is a more rigid statistical test of expression. This tool should
# be run on unnormalized data!)
# INPUT CDNA microarray[...].tsv OUTPUT singlechip.tsv
# PARAMETER do.lowess.normalization [yes, no] DEFAULT yes (Loess-normalize before analysis)
# PARAMETER analysis.method [noise-envelope, Newton] DEFAULT noise-envelope (Analysis method)
# PARAMETER significance.cutoff DECIMAL FROM 0 TO 10 DEFAULT 2 (Cut-off for significance)

# Single slide methods for cDNA chips
# JTT 4.7.2006

# Loading the libraries
library(sma)

# Renaming variables
norm<-do.lowess.normalization
meth<-analysis.method
cutoff<-significance.cutoff

# Loading the data
file<-dir()
dat<-read.table(file[1], header=T, sep="\t")

# Calculating A and M
x<-dat$sample-dat$samplebg
y<-dat$control+dat$controlbg
M<-log2(y/x)
A<-log2((x+y)/2)

# Normalization
if(norm=="no") {
   A<-A
   M<-M
}
if(norm=="yes") {
   ind <- is.na(A) | is.na(M) | is.infinite(A) | is.infinite(M)
   l <- lowess(A[!ind], M[!ind])
   ratio <- M
   ratio[!ind] <- M[!ind] - approx(l, xout = A[!ind])$y
   M<-ratio
}

# Analysis
if(meth=="noise-envelope") {
   b<-length(A)-25
   m<-c()
   sdm<-c()
   for(i in 1:b) {
      m<-c(m, mean(M[i:(i+25)]))
      sdm<-c(sdm, sd(M[i:(i+25)]))
   }
   m<-c(m, rep(tail(m, n=1), 25))
   sdm<-c(sdm, rep(tail(sdm, n=1), 25))
   res<-na.omit(ifelse(M>=(m+cutoff*sdm) | M<=(m-cutoff*sdm), M, NA))
   ind<-which(M>=(m+cutoff*sdm) | M<=(m-cutoff*sdm))
}

if(meth=="Newton") {
   # Newton function as implemented in sma package
   ind <- is.na(M) | is.na(A)
   theta <- func.em(A[!ind], M[!ind])
   logodds <- rep(NA, length(M))
   logodds[!ind] <- lod2(A[!ind], M[!ind], theta)
   res<-na.omit(ifelse(logodds>=cutoff | logodds<=-(cutoff), M, NA)) 
   ind<-which(logodds<=-cutoff | logodds>=cutoff)
}

# Writes a table of significant genes
write.table(data.frame(genes=dat$identifier[ind], chip.file1=M[ind], average.file1=A[ind], logodds=logodds[ind]), "singlechip.tsv", sep="\t", col.names=T, quote=F, row.names=F)
