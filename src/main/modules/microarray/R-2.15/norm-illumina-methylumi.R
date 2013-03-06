# TOOL norm-illumina-methylumi.R: "Illumina - methylumi pipeline" (Illumina methylation assay normalization using FinalReport files and lumi methodology. You need to import the FinalReport file DIRECTLY, not using the Import tool. The FinalReport file has to contain sample methylation profiles, group profile will not work with this tool.)
# INPUT FinalReport_sample_methylation_profile.txt: FinalReport_sample_methylation_profile.txt TYPE GENERIC 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT unmethylated.tsv: unmethylated.tsv 
# OUTPUT methylated.tsv: methylated.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# OUTPUT OPTIONAL QC-plot.pdf: QC-plot.pdf
# PARAMETER chiptype: Chiptype TYPE [HumanMethylation27: HumanMethylation27, HumanMethylation450: HumanMethylation450] DEFAULT HumanMethylation27 (Select the correct BeadChip type)
# PARAMETER OPTIONAL normalization: Normalization TYPE [none: none, quantile: quantile, ssn: ssn] DEFAULT quantile ()
# PARAMETER OPTIONAL color.balance.adjustment: "Color balance adjustment" TYPE [none: none, quantile: quantile, ssn: ssn] DEFAULT quantile (Adjustment of color balance)
# PARAMETER OPTIONAL background.correction: "Background correction" TYPE [none: none, bgAdjust2C: bgAdjust2C, forcePositive: forcePositive] DEFAULT none (Should background adjustment be applied)
# PARAMETER OPTIONAL QCplots: QCplots TYPE [yes: yes, no: no] DEFAULT yes (Do you want quality control plots)
# PARAMETER OPTIONAL image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the QC image)
# PARAMETER OPTIONAL image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the QC image)


# Illumina methylation array data preprocessing and normalization for FinalReport file
# JTT 2.2.2011
# Modified 25.9.2012 and 28.10.2012 by JTT

# setwd("C://Users//Jarno Tuimala//Desktop//methylumi data")
# color.balance.adjustment<-c("quantile")
# background.correction<-c("none")
# normalization<-c("quantile")
# chiptype<-c("HumanMethylation450")
# QCplots<-"yes"
# image.width<-c(600)
# image.height<-c(600)

# Renaming variables
w<-image.width
h<-image.height

# Loading libraries
library(lumi)
library(methylumi)
library(annotate)

# Converting to the correct chiptype
if(chiptype=="HumanMethylation27") {
	chiptype<-c("IlluminaHumanMethylation27k")
}
if(chiptype=="HumanMethylation450") {
	chiptype<-c("IlluminaHumanMethylation450k")
}
chiptype<-paste(chiptype, ".db", sep="")

# Loading data files
#dat<-lumiMethyR("chip.tsv", lib=chiptype)
dat<-methylumiR("FinalReport_sample_methylation_profile.txt", lib=chiptype)
methyLumiM <- as(dat, "MethyLumiM")
methyLumiM <- addAnnotationInfo(methyLumiM, lib = chiptype)

# Color balance adjustment
dat2<-lumiMethyC(methyLumiM, method=color.balance.adjustment)

# Background adjustment
if(color.balance.adjustment!="none") {
   dat3<-lumiMethyB(dat2, method=background.correction, separateColor=FALSE)
} else {
   dat3<-lumiMethyB(dat, method=background.correction, separateColor=TRUE)
}

# Normalization
dat4<-lumiMethyN(dat3, method=normalization)
#dat4<-normalizeMethyLumiSet(dat)

# QC plots
if(QCplots=="yes") {
   pdf(file="QC-plot.pdf", width=w/72, height=h/72)
   par(mfrow=c(2,2))
   plotColorBias1D(dat, main="Unpreprocessed")
   plotColorBias1D(dat4, main="Preprocessed")
   boxplotColorBias(dat, main="Unpreprocessed")
   boxplotColorBias(dat4, main="Preprocessed")
   dev.off()
}

# Convert sample names to Chipster style
dat5<-exprs(dat4)
sample.names<-colnames(dat5)
sample.names<-paste("chip.", sample.names, sep="")
names(dat5)<-sample.names
colnames(dat5)<-sample.names

dat6<-methylated(dat4)
sample.names<-colnames(dat6)
sample.names<-paste("chip.", sample.names, sep="")
names(dat6)<-sample.names
colnames(dat6)<-sample.names

dat7<-unmethylated(dat4)
sample.names<-colnames(dat7)
sample.names<-paste("chip.", sample.names, sep="")
names(dat7)<-sample.names
colnames(dat7)<-sample.names

# Annotations
library(chiptype, character.only=T)
symbols <- unlist (lookUp(rownames(dat5), chiptype, what="SYMBOL"))
genenames <- unlist (lookUp(rownames(dat5), chiptype, what="GENENAME"))
symbols <- gsub("#", "", symbols)
genenames <- gsub("#", "", genenames)
symbols <- gsub("'", "", symbols)
genenames <- gsub("'", "", genenames)

# Write out a phenodata
group<-c(rep("", ncol(dat5)))
training<-c(rep("", ncol(dat5)))
write.table(data.frame(sample=sample.names, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Write out expression data
write.table(data.frame(symbol=symbols, description=genenames, dat5), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
write.table(data.frame(symbol=symbols, description=genenames, dat6), file="methylated.tsv", col.names=T, quote=F, sep="\t", row.names=T)
write.table(data.frame(symbol=symbols, description=genenames, dat7), file="unmethylated.tsv", col.names=T, quote=F, sep="\t", row.names=T)
