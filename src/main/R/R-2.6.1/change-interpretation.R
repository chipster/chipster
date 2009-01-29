# ANALYSIS Utilities/"Change interpretation" (Let's user to transform the expression values.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT change-interpretation.tsv
# PARAMETER transform [log2-linear, linear-log2, log10-linear, linear-log10, ln-linear, linear-ln] DEFAULT log2-linear (From which to transform to what)


# Change interpretation
# JTT 26.4.2008

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat<-dat[,grep("chip", names(dat))]

# Transforms the data
if(transform=="log2-linear") {
   dat2<-2^dat
}
if(transform=="linear-log2") {
   if(any(dat<0)) {
      stop("Negative values in the data! Can't log-transform.")
   }
   dat2<-log2(dat)
}

if(transform=="log10-linear") {
   dat2<-10^dat
}
if(transform=="linear-log10") {
   if(any(dat<0)) {
      stop("Negative values in the data! Can't log-transform.")
   }
   dat2<-log10(dat)
}

if(transform=="ln-linear") {
   dat2<-exp(dat)
}
if(transform=="linear-ln") {
   if(any(dat<0)) {
      stop("Negative values in the data! Can't log-transform.")
   }
   dat2<-ln(dat)
}

# Writing the data to disk
write.table(dat2, "change-interpretation.tsv", sep="\t", row.names=T, col.names=T, quote=F)

