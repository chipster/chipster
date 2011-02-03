# ANALYSIS Normalisation/"Illumina - methylumi pipeline" (Illumina methylation assay normalization using FinalReport files and lumi methodology.
# TO USE THIS, YOU NEED TO IMPORT THE FinalReport FILE DIRECTLY, NOT USING THE IMPORT TOOL.)
# INPUT GENERIC chip.tsv OUTPUT normalized.tsv, unmethylated.tsv, methylated.tsv, phenodata.tsv
# PARAMETER color.balance.adjustment [none, quantile, ssn] DEFAULT quantile (Adjustment of color balance)
# PARAMETER background.correction [none, bgAdjust2C, forcePositive] DEFAULT none (Should background adjustment be applied)
# PARAMETER normalization [none, quantile, ssn] DEFAULT quantile (Normalization)
# PARAMETER chiptype [Human] DEFAULT Human (chiptype)


# Illumina methylation array data preprocessing and normalization for FinalReport file
# JTT 2.2.2011

# color.balance.adjustment<-c("quantile")
# background.correction<-c("none")
# normalization<-c("quantile")
# chiptype<-c("Human")


# Loading libraries
library(lumi)
library(methylumi)

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


# Write out a phenodata
group<-c(rep("", ncol(dat5)))
training<-c(rep("", ncol(dat5)))
write.table(data.frame(sample=sample.names, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

write.table(data.frame(dat5), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
write.table(data.frame(dat6), file="methylated.tsv", col.names=T, quote=F, sep="\t", row.names=T)
write.table(data.frame(dat7), file="unmethylated.tsv", col.names=T, quote=F, sep="\t", row.names=T)
