# TOOL norm-agilent-miRNA.R: "Agilent miRNA" (Agilent miRNA chip data preprocessing. Automatically averages all the rows,
# i.e., miRNA:s that have the same name. YOU HAVE TO SPECIFY THE CHIPTYPE.
# To be able to remove the control probes, the column labelled "ControlType" in the raw data file should be marked as
# "Annotation" during importing.)
# INPUT microarray{…}.tsv: "Raw data files" TYPE CDNA
# OUTPUT normalized.tsv: "Normalized data"
# OUTPUT phenodata.tsv: "Experiment description"
# PARAMETER background.treatment: "Background treatment" TYPE [none, subtract, edwards, normexp] DEFAULT normexp (Background treatment method)
# PARAMETER background.offset: "Background offset" TYPE [0, 50] DEFAULT 50 (Background offset)
# PARAMETER normalize.chips: "Normalize chips" TYPE [none, scale, scale-75, quantile, vsn] DEFAULT quantile (Between arrays normalization method)
# PARAMETER remove.control.probes: "Remove control probes" TYPE [yes, no] DEFAULT no (Remove control probes from the dataset)
# PARAMETER chiptype: "Chiptype" TYPE [empty, "Human", "Mouse", "Rat"]  DEFAULT empty (Chiptype)

# Agilent miRNA chip normalization
# MG 
# 21.10.2010


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
if(normba=="scale-75") {
   dat3<-normalizeBetweenArrays(dat2$R, method="none")
   normfact<-apply(dat3, 2, quantile)["75%",]
   for(i in 1:ncol(dat3)) {
      dat3[,i]<-dat3[,i]/normfact[i]
   }
} else {
   dat3<-normalizeBetweenArrays(dat2$R, method=normba)
}

# Log-transforming the data
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

# Annotation package name conversions}
if(chiptype=="Human") {
	chiptype<-c("miRNA")
}
if(chiptype=="Mouse") {
	chiptype<-c("miRNA")
}
if(chiptype=="Rat") {
	chiptype<-c("miRNA")
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
names(M)<-paste(names(M), ".tsv", sep="")
M<-data.frame(M)
rownames(M)<-genes

if(chiptype=="cDNA") {
   if(nrow(flags)!=nrow(M)) {
      write.table(data.frame(round(M, digits=2)), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
   } else {
     write.table(data.frame(round(M, digits=2), flags), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
   }
}

if(chiptype=="miRNA") {
	if(nrow(flags)!=nrow(M)) {
		write.table(data.frame(round(M, digits=2)), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	} else {
		write.table(data.frame(round(M, digits=2), flags), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	}
}

