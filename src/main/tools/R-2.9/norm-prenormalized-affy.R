# ANALYSIS Normalisation/"Process prenormalized affy" (If you import prenormalized Affymetrix data that is not in Chipster format, you
# need to import it through Import tool, and then use this tool for preprocessing it. Every column containing normalized
# expression values need to be specifies as Sample in the Import tool. If you want to be able to use annotation,
# you need to SPECIFY THE CHIPTYPE.)
# INPUT CDNA microarray[...].tsv OUTPUT normalized.tsv, phenodata.tsv
# PARAMETER chiptype STRING DEFAULT empty (chiptype)

# Process prenormalized affy
# MG 21.10.2009

# Loads the libraries
library(limma)
library(affy)
library(gcrma)

# Reading data

columns<-list(R="sample", Rb="sample", G="sample", Gb="sample")
annotation<-c("identifier")
columns.other<-c("flag", "annotation")

files<-dir()
files<-files[files!="phenodata.tsv"]
dat<-read.maimages(files=files, columns=columns, annotation=annotation, other.columns=columns.other) 

#dat <- read.table(file="normalized.tsvt", header=TRUE, sep="\t")

# Extract annotations
#if (keep.annotations=="yes") {
#	annotation.info <- dat$other[[1]]
#	annotation.info <- as.character (annotation.info[,1])
#}

# Mock normalization
dat2<-normalizeBetweenArrays(dat$R, method="none")

# Writes out a phenodata table
sample<-paste(colnames(dat2), ".tsv", sep="")
group<-c(rep("", length(sample)))
training<-c(rep("", length(sample)))
time<-c(rep("", length(sample)))
random<-c(rep("", length(sample)))
if(chiptype=="empty") {
	chiptype<-c("cDNA")
}
write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)


# Preparing data for export
M<-data.frame(dat2)
colnames(M)<-paste("chip.", colnames(M), sep="")
M2<-aggregate(M, as.list(dat$genes), mean)
rownames(M2)<-M2$identifier
M2<-M2[,-1]

# Writing out data
a<-try(library(chiptype, character.only=T))
if(chiptype!="empty" & class(a)!="try-error") {
	# Including gene names to data
	lib2<-sub('.db','',chiptype)
	symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep="")))))[rownames(M2),])
	genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep="")))))[rownames(M2),])
	# Fxes an issue introduced in BioC2.4 where the "#" character is introduced in some gene names
	genename <- gsub("#", "", genename)
	# Writes the results into a file
	write.table(data.frame(symbol, description=genename, M2, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T))
} 

if(chiptype=="empty" | class(a)=="try-error") {
	write.table(data.frame(M2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
}

