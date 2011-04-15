# ANALYSIS Statistics/"Time series" (Analyses of time series data. Finds periodically expressed genes. 
# For the ICA method, a standard deviation threshold for significantly differentially expressed genes is needed.
# If there are more than one replicate per time point, this tool will not work.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT timeseries.tsv, profiles.pdf
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the time to test)
# PARAMETER analysis.type [periodicity, ica] DEFAULT periodicity (Analysis type)
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER p.value.adjustment.method [yes, no] DEFAULT yes (Apply Benjamimi-Hochberg correction?)
# PARAMETER SD.for.ICA DECIMAL FROM 0 TO 10 DEFAULT 2.0 (Standard deviation for ICA)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Analysis methods for timeseries
# JTT 21.7.2006

# Parameter settings (default) for testing purposes
#column<-"time"
#analysis.type<-"periodicity"
#p.value.threshold<-0.05
#p.value.adjustment.method<-"yes"
#SD.for.ICA<-2.0
#image.width<-600
#image.height<-600

# Load the libraries
library(e1071)
library(GeneCycle)
library(fastICA)

# Renaming variables
analysis<-analysis.type
p.cut<-p.value.threshold
multiple.correction<-p.value.adjustment.method
thresh<-SD.for.ICA
w<-image.width 
h<-image.height 

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Reads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
times<-phenodata[,grep(column, colnames(phenodata))]

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# How many replicates there are per time point?
if(length(times)>length(levels(as.factor(times)))) {
   repl<-c()
   for(i in 1:length(levels(as.factor(times)))) {
      repl<-c(repl, sum(as.numeric(grep(times[1], times, value=T)==times[1])))
   }
} else {
   repl<-rep(1, length(times))
}

# Making a longitudinal object
dat3<-as.longitudinal(t(dat2), repeats=repl, time=times)

# Replacing missing values
dat4<-t(na.omit(t(impute(dat3))))

if(analysis=="periodicity") {
   f<-fisher.g.test(dat4)
   if(multiple.correction=="yes") {
      # Estimates the proportion of null p-values, uses BH-method
      p.adj<-fdrtool(f)$lfdr
      dev.off()
   } else {
      p.adj<-f
   }
   dat5<-as.data.frame(t(dat4))
   names(dat5)<-names(dat2)
   write.table(data.frame(dat5[p.adj<=p.cut,], p.adjusted=round(p.adj[p.adj<=p.cut], digits=6)), file="timeseries.tsv", sep="\t", row.names=T, col.names=T, quote=F)
   pdf(file="profiles.pdf", width=w/72, height=h/72)
   plot(1, 1, col=0)
   text(1, 1, "This is a dummy image.", col=1)
   text(1, 0.9, "To generate an image of gene expression profiles, use ica option.", col=1)
   dev.off()
}


if(analysis=="ica") {
   # Calculating independent component analysis usign a fast method
   o<-fastICA(t(dat4), n.comp=(ncol(t(dat4))-1), method="C")
   d<-c()
   g<-c()
   for(i in 1:ncol(o$S)) {
      s<-o$S[,i]
      l<-which(s>=(mean(s)+thresh*sd(s)))
      d<-c(d, l)
      g<-c(g, rep(i, length(l)))
      dg<-data.frame(dat[d,], cluster=g)
   }
   write.table(dg, file="timeseries.tsv", sep="\t", row.names=T, col.names=T, quote=F)
   a<-data.frame(times, t(o$A))
   pdf(file="profiles.pdf", width=w/72, height=h/72)
   par(mar=c(0, 1, 1, 0)+0.1)
   par(mfrow=c(ceiling(sqrt(ncol(a))), ceiling(sqrt(ncol(a)))))
   for(i in 2:ncol(a)) {
      plot(a$time, a[,i], xaxt="n", type="l", xlab=NULL, ylab=NULL, main=paste("Chip", i-1, sep=" "), cex.main=0.75, cex.axis=0.7, cex.lab=0.75, tck=-0.01, mgp=c(3,0.2,0))
   }
   dev.off()
}
