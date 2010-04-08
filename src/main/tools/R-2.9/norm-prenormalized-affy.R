# ANALYSIS Normalisation/"Process prenormalized affy" (If you import prenormalized Affymetrix data that is not in Chipster format, you
# need to import it through the Import Wizard, and then use this tool for preprocessing it. During data import, make sure to mark every column containing normalized
# expression values as Sample and the column containing the Affymetrix probe ID:s as Identifier. If you want to be able to use annotation,
# you need to SPECIFY THE CHIPTYPE, e.g. hgu133a2.db.)
# INPUT CDNA microarray[...].tsv OUTPUT normalized.tsv, phenodata.tsv
# PARAMETER chiptype STRING DEFAULT empty (chiptype)

# Process prenormalized affy
# MG 6.4.2010

# Loads the libraries
library(limma)
library(affy)
library(gcrma)

# Reading data

columns<-list(R="sample")
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
sample<-colnames(dat2)
colnames(dat2)<-paste("chip.", colnames(dat2), sep="")

# Writes out a phenodata table
# sample<-paste(sample, ".tsv", sep="")
group<-c(rep("", length(sample)))
training<-c(rep("", length(sample)))
time<-c(rep("", length(sample)))
random<-c(rep("", length(sample)))
if(chiptype=="empty") {
	chiptype<-"Affy"
}
write.table(data.frame(sample=sample, chiptype=chiptype, group=group, description=sample), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)


# Writing out data
a<-try(library(chiptype, character.only=T))
if(chiptype!="empty" & class(a)!="try-error") {
	# Including gene names to data
	lib2<-sub('.db','',chiptype)
	symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep="")))))[unlist(dat$genes),])
	genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep="")))))[unlist(dat$genes),])
	# Fxes an issue introduced in BioC2.4 where the "#" character is introduced in some gene names
	genename <- gsub("#", "", genename)
	# Writes the results into a file
	output_table <- data.frame(symbol, description=genename, dat2)
	rownames (output_table) <- unlist(dat$genes)
	write.table(output_table, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
} 

if(chiptype=="empty" | class(a)=="try-error") {
	output_table <- data.frame(dat2)
	rownames (output_table) <- unlist(dat$genes)
	write.table(output_table, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
}

