# TOOL norm-illumina-methylumi.R: "Illumina - methylumi pipeline" (Illumina methylation assay normalization using FinalReport files and lumi methodology. TO USE THIS, YOU NEED TO IMPORT THE FinalReport FILE DIRECTLY, NOT USING THE IMPORT TOOL.)
# INPUT chip.tsv: chip.tsv TYPE GENERIC 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT unmethylated.tsv: unmethylated.tsv 
# OUTPUT methylated.tsv: methylated.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# OUTPUT QC-plot.pdf: QC-plot.pdf 
# PARAMETER color.balance.adjustment: color.balance.adjustment TYPE [none: none, quantile: quantile, ssn: ssn] DEFAULT quantile (Adjustment of color balance)
# PARAMETER background.correction: background.correction TYPE [none: none, bgAdjust2C: bgAdjust2C, forcePositive: forcePositive] DEFAULT none (Should background adjustment be applied)
# PARAMETER normalization: normalization TYPE [none: none, quantile: quantile, ssn: ssn] DEFAULT quantile ()
# PARAMETER chiptype: chiptype TYPE [Human: Human] DEFAULT Human ()
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the QC image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the QC image)


# Illumina methylation array data preprocessing and normalization for FinalReport file
# JTT 2.2.2011

# color.balance.adjustment<-c("quantile")
# background.correction<-c("none")
# normalization<-c("quantile")
# chiptype<-c("Human")
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
if(chiptype=="Human") {
	chiptype<-c("IlluminaHumanMethylation27k")
}
chiptype<-paste(chiptype, ".db", sep="")

# Loading data files
dat<-lumiMethyR("chip.tsv", lib="IlluminaHumanMethylation27k.db")

# Color balance adjustment
dat2<-lumiMethyC(dat, method=color.balance.adjustment)

# Background adjustment
if(color.balance.adjustment!="none") {
   dat3<-lumiMethyB(dat2, method=background.correction, separateColor=FALSE)
} else {
   dat3<-lumiMethyB(dat, method=background.correction, separateColor=TRUE)
}

# Normalization
dat4<-lumiMethyN(dat3, method=normalization)

# QC plots
pdf(file="QC-plot.pdf", width=w/72, height=h/72)
par(mfrow=c(2,2))
plotColorBias1D(dat, main="Unpreprocessed")
plotColorBias1D(dat4, main="Preprocessed")
boxplotColorBias(dat, main="Unpreprocessed")
boxplotColorBias(dat4, main="Preprocessed")
dev.off()

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
