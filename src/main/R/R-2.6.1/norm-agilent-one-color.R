# ANALYSIS Normalisation/"Agilent 1-color" (Agilent one-color data preprocessing. Automatically averages all the rows,
# i.e., genes that have the same name. YOU HAVE TO SPECIFY THE CHIPTYPE.)
# INPUT CDNA microarray[...].tsv OUTPUT normalized.tsv, phenodata.tsv
# PARAMETER background.treatment [none, subtract, edwards, normexp] DEFAULT normexp (Background treatment method)
# PARAMETER background.offset [0, 50] DEFAULT 50 (Background offset)
# PARAMETER normalize.chips [none, scale, quantile, vsn] DEFAULT none (Between arrays normalization method)
# PARAMETER remove.control.probes [yes, no] DEFAULT no (Remove control probes from the dataset)
# PARAMETER chiptype [empty, Human-1 (4100a), Human-2 (4101a), Human-1A (4110b), Human-1B (4111a), Human-Whole-Genome (4112a), Mouse (4104a), Mouse (4120a), Mouse (4121a), Mouse (4122a), Rat (4105a), Rat (4130a), Rat (4131), zebrafish, zebrafishV2, drosophila, rhesus, rice] DEFAULT empty (chiptype)


# cDNA chip normalization
# JTT 15.10.2007

# Loads the libraries
library(limma)

# Renaming variables
bg<-background.treatment
normba<-normalize.chips

# Reading data
columns<-list(R="sample", Rb="samplebg", G="sample", Gb="samplebg")
annotation<-c("identifier")
columns.other<-c("flag", "annotation")

files<-dir()
files<-files[files!="phenodata.tsv"]
dat<-read.maimages(files=files, columns=columns, annotation=annotation, other.columns=columns.other) 

# Background correction
dat2<-backgroundCorrect(dat, bg, offset=as.numeric(background.offset))

# Normalization across arrays
dat3<-normalizeBetweenArrays(dat2$R, method=normba)
dat3<-log2(dat3)

# Writes out a phenodata table
sample<-paste(colnames(dat3), ".tsv", sep="")
group<-c(rep("", length(sample)))
training<-c(rep("", length(sample)))
time<-c(rep("", length(sample)))
random<-c(rep("", length(sample)))
if(chiptype=="empty") {
   chiptype<-c("cDNA")
}
if(chiptype=="Human-1(4100a)") {
   chiptype<-c("hgug4100a")
}
if(chiptype=="Human-2(4101a)") {
   chiptype<-c("hgug4101a")
}
if(chiptype=="Human-1A(4110b)") {
   chiptype<-c("hgug4110b")
}
if(chiptype=="Human-1B(4111a)") {
   chiptype<-c("hgug4111a")
}
if(chiptype=="Human-Whole-Genome(4112a)") {
   chiptype<-c("hgug4112a")
}
if(chiptype=="Mouse(4104a)") {
   chiptype<-c("mgug4104a")
}
if(chiptype=="Mouse(4120a)") {
   chiptype<-c("mgug4120a")
}
if(chiptype=="Mouse(4121a)") {
   chiptype<-c("mgug4121a")
}
if(chiptype=="Mouse(4122a)") {
   chiptype<-c("mgug4122a")
}
if(chiptype=="Rat(4105a)") {
   chiptype<-c("rgug4105a")
}
if(chiptype=="Rat(4130a)") {
   chiptype<-c("rgug4130a")
}
if(chiptype=="Rat(4131)") {
   chiptype<-c("rgug4131unigene")
}
if(chiptype=="zebrafish") {
   chiptype<-c("AgilentZF")
}
if(chiptype=="zebrafishV2") {
   chiptype<-c("AZFv2")
}
if(chiptype=="drosophila") {
   chiptype<-c("DrosoAgilent")
}
if(chiptype=="rhesus") {
   chiptype<-c("rhesusrefseq")
}
if(chiptype=="rice") {
   chiptype<-c("rice")
}
write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Removes control probes
rownames(dat3)<-dat$genes$identifier
if(remove.control.probes=="yes") {
   if(is.null(dim(dat$other$annotation))==FALSE) {
      dat3<-dat3[rowSums(dat$other$annotation)==0,]
   }
}

# Constructs and writes out a table
M<-dat3
M<-aggregate(M, list(rownames(M)), mean)
rownames(M)<-M$Group.1
genes<-rownames(M)
M<-M[,-1]
if(length(dat$other$flag)!=0) {
   flags<-as.data.frame(dat$other$flag)
   names(flags)<-paste("flag.", names(flags), sep="")
}
if(length(dat$other$flag)==0) {
   flags<-matrix(nrow=0, ncol=0)
}
names(M)<-paste("chip.", names(M), sep="")
M<-data.frame(M)
rownames(M)<-genes

if(chiptype!="cDNA") {
   # Including gene names to data
   library(chiptype, character.only=T)
   symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(chiptype, "SYMBOL", sep="")))))[rownames(M),])
   genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(chiptype, "GENENAME", sep="")))))[rownames(M),])
}

if(chiptype!="cDNA") {
   # Conditional on whether flags are included or not, write the data to disk
   if(nrow(flags)!=nrow(M)) {
      write.table(data.frame(symbol, description=genename, round(M, digits=2)), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
   } else {
     write.table(data.frame(symbol, description=genename, round(M, digits=2), flags), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
   }
}

if(chiptype=="cDNA") {
   if(nrow(flags)!=nrow(M)) {
      write.table(data.frame(round(M, digits=2)), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
   } else {
     write.table(data.frame(round(M, digits=2), flags), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
   }
}
