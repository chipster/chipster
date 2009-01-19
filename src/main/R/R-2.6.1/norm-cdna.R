# ANALYSIS Normalisation/cDNA (cDNA data preprocessing. Automatically averages all the rows,
# i.e., genes that have the same name.)
# INPUT CDNA microarray[...].tsv OUTPUT normalized.tsv, phenodata.tsv
# PARAMETER background.treatment [none, subtract, edwards, normexp] DEFAULT normexp (Background treatment method)
# PARAMETER background.offset [0, 50] DEFAULT 50 (Background offset)
# PARAMETER normalize.arrays [none, median, loess] DEFAULT loess (Within array normalization method)
# PARAMETER normalize.genes [none, scale, quantile, Aquantile, vsn] DEFAULT none (Between arrays normalization method)


# cDNA chip normalization
# JTT 9.6.2006

# Loads the libraries
library(limma)

# Renaming variables
bg<-background.treatment
normwa<-normalize.arrays
normba<-normalize.genes

# Reading data
columns<-list(R="sample", Rb="samplebg", G="control", Gb="controlbg")
annotation<-c("annotation", "identifier")
columns.other<-c("flag")

files<-dir()
files<-files[files!="phenodata.tsv"]
dat<-read.maimages(files=files, columns=columns, annotation=annotation, other.columns=columns.other) 

# Normalization within arrays
dat2<-normalizeWithinArrays(dat, method=normwa, bc.method=bg, offset=as.numeric(background.offset))

# Normalization across arrays
dat3<-normalizeBetweenArrays(dat2, method=normba)

# Modifies the phenodata table by writing one more row to it ((Affymetrix) chip type)
sample<-paste(colnames(dat2$M), ".tsv", sep="")
group<-c(rep("", length(sample)))
training<-c(rep("", length(sample)))
time<-c(rep("", length(sample)))
random<-c(rep("", length(sample)))
chiptype<-c("cDNA")
write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Constructs and writes out a table of data
M<-dat3$M
A<-dat3$A
rownames(M)<-dat3$genes$identifier
rownames(A)<-dat3$genes$identifier
M<-aggregate(M, list(rownames(M)), mean)
A<-aggregate(A, list(rownames(A)), mean)
rownames(M)<-M$Group.1
rownames(A)<-A$Group.1
genes<-rownames(M)
M<-M[,-1]
A<-A[,-1]
if(length(dat$other$flag)!=0) {
   flags<-as.data.frame(dat$other$flag)
   names(flags)<-paste("flag.", names(flags), sep="")
}
if(length(dat$other$flag)==0) {
   flags<-matrix(nrow=0, ncol=0)
}
names(M)<-paste("chip.", names(M), sep="")
names(A)<-paste("average.", names(A), sep="")

A<-data.frame(round(A, digits=2))
M<-data.frame(round(M, digits=2))

if(nrow(flags)!=nrow(A) | nrow(flags)!=nrow(M)) {
   D<-data.frame(M, A)
} else {
   D<-data.frame(M, A, flags)
}
rownames(D)<-genes

write.table(D, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
