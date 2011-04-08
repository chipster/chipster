# ANALYSIS Normalisation/"Illumina - lumi pipeline" (Illumina normalization using BeadSummaryData files and lumi methodology.
# TO USE THIS, YOU NEED TO IMPORT THE BeadSummaryData FILE DIRECTLY, NOT USING THE IMPORT TOOL.)
# INPUT GENERIC chip.tsv OUTPUT normalized.tsv, phenodata.tsv
# PARAMETER transformation [none, vst, log2] DEFAULT log2 (How to transform the data)
# PARAMETER background.correction [none, bgAdjust.affy] DEFAULT none (Should background adjustment be applied)
# PARAMETER normalize.chips [none, rsn, loess, quantile, vsn] DEFAULT quantile (Between arrays normalization method)
# PARAMETER chiptype [empty, Human, Mouse, Rat] DEFAULT empty (chiptype)


# Illumina data preprocessing and normalization for BeadSummaryData
# JTT 27.5.2008
#
# NG, 21.12.2010 modified to add gene symbol and gene name to the output

# transformation<-"none"
# normalize.genes<-"none"
# chiptype<-"Mouse-6v1.1"

# Loading libraries
library(lumi)
library(annotate)

# Converting to the correct chiptype
if(chiptype=="empty") {
	chiptype<-c("Illumina")
	mapping<-c("Illumina")
}
if(chiptype=="Human") {
	chiptype<-c("lumiHumanAll")
	mapping<-c("lumiHumanIDMapping")
}
if(chiptype=="Mouse") {
	chiptype<-c("lumiMouseAll")
	mapping<-c("lumiMouseIDMapping")
}
if(chiptype=="Rat") {
	chiptype<-c("lumiRatAll")
	mapping<-c("lumiRatIDMapping")
}
chiptype<-paste(chiptype, ".db", sep="")

# Loading data files
# fileName <-dir()
# fileName<-fileName[fileName!="phenodata.tsv"]
x.lumi <- lumiR("chip.tsv", lib.mapping=mapping)

# Quality control (not run)
# QC
# q.lumi <- lumiQ(x.lumi)  

# Background correction
if(background.correction=="bgAdjust.affy") {
   x.lumi<-lumiB(x.lumi, method="bgAdjust.affy")
}

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
if(normalize.chips=="none") {
   n.lumi<-t.lumi
}
if(normalize.chips=="rsn") {
   n.lumi<-lumiN(t.lumi, method=c("rsn"))
}
if(normalize.chips=="loess") {
   n.lumi<-lumiN(t.lumi, method=c("loess"))
}
if(normalize.chips=="quantile") {
   n.lumi<-lumiN(t.lumi, method=c("quantile"))
}
if(normalize.chips=="vsn") {
   n.lumi<-lumiN(t.lumi, method=c("vsn"))
}

# Convert sample names to Chipster style
dat2<-exprs(n.lumi)
sample.names<-colnames(dat2)
sample.names<-paste("chip.", sample.names, sep="")
names(dat2)<-sample.names
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
   # symbol<-gsub("#", "", symbol)
   # genename<-gsub("#", "", genename)
   # Write out expression data
   # write.table(data.frame(symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	symbols <- unlist (lookUp(rownames(dat2), chiptype, what="SYMBOL"))
	genenames <- unlist (lookUp(rownames(dat2), chiptype, what="GENENAME"))
	symbols <- gsub("#", "", symbols)
	genenames <- gsub("#", "", genenames)
	write.table(data.frame(symbol=symbols, description=genenames, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
} else {
   # Write out expression data
   write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
}