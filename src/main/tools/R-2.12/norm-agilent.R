# ANALYSIS Normalisation/"Agilent 2-color" (Agilent two-color data preprocessing. Automatically averages all the rows,
# i.e., genes that have the same name. YOU HAVE TO SPECIFY THE CHIPTYPE.)
# INPUT CDNA microarray[...].tsv OUTPUT normalized.tsv, phenodata.tsv
# PARAMETER background.treatment [none, subtract, edwards, normexp] DEFAULT normexp (Background treatment method)
# PARAMETER background.offset [0, 50] DEFAULT 50 (Background offset)
# PARAMETER normalize.arrays [none, median, loess] DEFAULT loess (Within array normalization method)
# PARAMETER normalize.genes [none, scale, quantile, Aquantile, vsn] DEFAULT none (Between arrays normalization method)
# PARAMETER remove.control.probes [yes, no] DEFAULT no (Remove control probes from the dataset)
# PARAMETER chiptype [empty, Human-1 (4100a), Human-2 (4101a), Human-1A (4110b), Human-1B (4111a), Human-Whole-Genome (4112a), Mouse (4104a), Mouse (4120a), Mouse (4121a), Mouse (4122a), Rat (4105a), Rat (4130a), Rat (4131)] DEFAULT empty (chiptype)

# cDNA chip normalization
# JTT 9.6.2006
#
# MG, 10.3.2010
# modified script to cope with manually set flags

#background.treatment<-"normexp"
#background.offset<-50
#normalize.arrays<-"loess"
#normalize.genes<-"none"
#remove.control.probes<-"no"
#chiptype<-"4101a"

# Loads the libraries
library(limma)

# Renaming variables
bg<-background.treatment
normwa<-normalize.arrays
normba<-normalize.genes

# Reading data
columns<-list(R="sample", Rb="samplebg", G="control", Gb="controlbg")
annotation<-c("identifier")
columns.other<-c("flag", "annotation")

files<-dir()
files<-files[files!="phenodata.tsv"]
dat<-read.maimages(files=files, columns=columns, annotation=annotation, other.columns=columns.other) 

# Normalization within arrays
dat2<-normalizeWithinArrays(dat, method=normwa, bc.method=bg, offset=as.numeric(background.offset))

# Normalization across arrays
dat3<-normalizeBetweenArrays(dat2, method=normba)

# Writes out a phenodata table
sample<-paste(colnames(dat2$M), ".tsv", sep="")
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
if(chiptype!="cDNA") {
	chiptype<-paste(chiptype, ".db", sep="")
}
# chiptype<-paste(chiptype, ".db", sep="")

write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Removes control probes
if(remove.control.probes=="yes") {
	if(is.null(dim(dat$other$annotation))==FALSE) {
		dat3<-dat3[rowSums(dat$other$annotation)==0,]
	}
}


# Define function for averaging flags
average_flags <- function (flags.vector) {
	number_probes <- length(flags.vector)
	if (max(match(flags.vector,"P",nomatch=0))>0) {
		return("P")
	} else {
		if (max(match(flags.vector,"M",nomatch=0))>0) {
			return("M")
		} else {
			return("A")
		}
	}
}			



# Constructs and writes out a table
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


# If flags are available aggregate those as well
if(length(dat$other$flag)!=0) {
	flags <- dat$other$flag
	if(remove.control.probes=="yes") {
		if(is.null(dim(dat$other$annotation))==FALSE) {
			flags <- flags[rowSums(dat$other$annotation)==0,]
		}
	}
	rownames(flags) <- dat3$genes$identifier
	flags2 <- aggregate(flags, list(rownames(flags)), average_flags)
	rownames(flags2) <- flags2$Group.1
	flags3 <- flags2 [,-1]
	flags4 <- as.data.frame(flags3)
	names(flags4)<-paste("flag.", names(flags4), sep="")
}
if(length(dat$other$flag)==0) {
	flags<-matrix(nrow=0, ncol=0)
}

names(M)<-paste("chip.", names(M), sep="")
names(M)<-paste(names(M), ".tsv", sep="")
names(A)<-paste("average.", names(A), sep="")
names(A)<-paste(names(A), ".tsv", sep="")
A<-data.frame(A)
M<-data.frame(M)
rownames(M)<-genes
rownames(A)<-genes

if(chiptype!="cDNA") {
	# Including gene names to data
	library(chiptype, character.only=T)
	lib2<-sub('.db','',chiptype)
	symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep="")))))[rownames(M),])
	genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep="")))))[rownames(M),])
	symbol<-gsub("#", "", symbol)
	genename<-gsub("#", "", genename)
	symbol <- gsub("'", "", symbol)
	genename <- gsub("'", "", genename)
}

if(chiptype!="cDNA") {
	# Conditional on whether flags are included or not, write the data to disk
	if(length(dat$other$flag)==0) {
		write.table(data.frame(symbol, description=genename, round(M, digits=2), round(A, digits=2)), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	} else {
		write.table(data.frame(symbol, description=genename, round(M, digits=2), round(A, digits=2), flags4), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	}
}

if(chiptype=="cDNA") {
	if(length(dat$other$flag)==0) {
		write.table(data.frame(round(M, digits=2), round(A, digits=2)), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	} else {
		write.table(data.frame(round(M, digits=2), round(A, digits=2), flags4), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	}
}
