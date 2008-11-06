# ANALYSIS Normalisation/"Illumina - lumi pipeline" (Illumina normalization using BeadSummaryData files and lumi methodology.
# TO USE THIS, YOU NEED TO IMPORT THE BeadSummaryData FILE DIRECTLY, NOT USING THE IMPORT TOOL.)
# INPUT GENERIC chip.tsv OUTPUT normalized.tsv, phenodata.tsv
# PARAMETER transformation [none, vst, log2] DEFAULT log2 (How to transform the data)
# PARAMETER normalize.genes [none, rsn, loess, quantile, vsn] DEFAULT quantile (Between arrays normalization method)
# PARAMETER chiptype [empty, Human-6v1, HumanRef-8v1, Human-6v2, HumanRef-8v2, Mouse-6v1.0a, MouseRef-8v1.0a, RatRef-12] DEFAULT empty (chiptype)


# Illumina data preprocessing and normalization for BeadSummaryData
# JTT 27.5.2008

# transformation<-"none"
# normalize.genes<-"none"
# chiptype<-"Mouse-6v1.1"

# Loading libraries
library(lumi)

# Converting to the correct chiptype
   if(chiptype=="empty") {
      chiptype<-c("Illumina")
   }
   if(chiptype=="Human-6v1" | chiptype=="HumanRef-8v1") {
      chiptype<-c("lumiHumanV1")
   }
   if(chiptype=="Human-6v2" | chiptype=="HumanRef-8v2") {
      chiptype<-c("lumiHumanV2")
   }
   if(chiptype=="Mouse-6v1.0a" | chiptype=="MouseRef-8v1.0a") {
      chiptype<-c("lumiMouseV1")
   }
   if(chiptype=="RatRef-12") {
      chiptype<-c("lumiRatV1")
   }

# Loading data files
# fileName <-dir()
# fileName<-fileName[fileName!="phenodata.tsv"]
x.lumi <- lumiR("chip.tsv", lib=chiptype)

# Quality control (not run)
# QC
# q.lumi <- lumiQ(x.lumi)  

# Transformation
if(transformation=="none") {
   t.lumi<-x.lumi
}
if(transformation=="vst") {
   t.lumi<-lumiT(x.lumi, method=c("vst"))
}
if(transformation=="log2") {
   t.lumi<-lumiT(x.lumi, method=c("log2"))
}

# Normalization
if(normalize.genes=="none") {
   n.lumi<-t.lumi
}
if(normalize.genes=="rsn") {
   n.lumi<-lumiN(t.lumi, method=c("rsn"))
}
if(normalize.genes=="loess") {
   n.lumi<-lumiN(t.lumi, method=c("loess"))
}
if(normalize.genes=="quantile") {
   n.lumi<-lumiN(t.lumi, method=c("quantile"))
}
if(normalize.genes=="vsn") {
   n.lumi<-lumiN(t.lumi, method=c("vsn"))
}

# Convert sample names to Chipster style
dat2<-exprs(n.lumi)
sample.names<-colnames(dat2)
sample.names<-paste("chip.", sample.names, sep="")
names(dat2)<sample.names
colnames(dat2)<-sample.names

# Write out a phenodata
group<-c(rep("", ncol(dat2)))
training<-c(rep("", ncol(dat2)))
write.table(data.frame(sample=sample.names, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

if(chiptype!="Illumina") {
   # Including gene names to data
   library(chiptype, character.only=T)
   # symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(chiptype, "SYMBOL", sep="")))))[x.lumi@featureData@data$TargetID,])
   # genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(chiptype, "GENENAME", sep="")))))[x.lumi@featureData@data$TargetID,])
   # Write out expression data
   # write.table(data.frame(symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
   write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
} else {
   # Write out expression data
   write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
}