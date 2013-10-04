# TOOL norm-prenormalized.R: "Process prenormalized" (If you import prenormalized data that is not in Chipster format, you need to import it through Import tool, and then use this tool for preprocessing it. Every column containing normalized expression values need to be specified as Sample in the Import tool. If you want to be able to use annotation, you need to SPECIFY THE CHIPTYPE.)
# INPUT microarray{...}.tsv: microarray{...}.tsv TYPE CDNA 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER chiptype: chiptype TYPE STRING DEFAULT empty ()
# PARAMETER keep.annotations: keep.annotations TYPE [yes: yes, no: no] DEFAULT no (Keep or discard annotation column after preprocessing)


# Process prenormalized
# JTT 26.01.2009
# MG 21.10.2009
# EK 23.04.2013, added .tsv ending to expression column names so that sample renaming is possible in interactive visualizations
# MK 04.10.2013, special characters trimmed off

# Loads the libraries
library(limma)

# Reading data
columns<-list(R="sample", Rb="sample", G="sample", Gb="sample")
annotation<-c("identifier")
columns.other<-c("flag", "annotation")

files<-dir()
files<-files[files!="phenodata.tsv"]
dat<-read.maimages(files=files, columns=columns, annotation=annotation, other.columns=columns.other) 

# Extract annotations
if (keep.annotations=="yes") {
	annotation.info <- dat$other[[1]]
	annotation.info <- as.character (annotation.info[,1])
}

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
colnames(M)<-paste("chip.", colnames(M), ".tsv", sep="")
M2<-aggregate(M, as.list(dat$genes), mean)
rownames(M2)<-M2$identifier
M2<-M2[,-1]

# Match annotations to identifiers

if (keep.annotations=="yes") {
	matching_table <- cbind (as.character (dat$genes[,1]), annotation.info)
	rownames (matching_table) <- as.character (dat$genes[,1])
	matching_table <- as.data.frame (matching_table)
	annotations <- as.character (matching_table[rownames(M2),2])
}

if(chiptype!="cDNA") {
	# Including gene names to data
	lib2<-sub('.db','',chiptype)
	library(chiptype, character.only=T)
	symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep="")))))[rownames(M2),])
	genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep="")))))[rownames(M2),])	
	genename <- gsub("#", "", genename)
	dat2 <- data.frame(symbol, description=genename, round(M2, digits=2))
}

if(chiptype=="cDNA" & keep.annotations=="no") {
	dat2 < data.frame(round(M2, digits=8))
}
if(chiptype=="cDNA" & keep.annotations=="yes") {
	dat2 <- data.frame(annotations, round(M2, digits=8))
}

if(length(grep("description", tolower(colnames(dat2)))) > 0) {
	dat2[, grep("description", tolower(colnames(dat2)))] <- gsub("\'+", "", dat2[, grep("description", tolower(colnames(dat2)))])
	dat2[, grep("description", tolower(colnames(dat2)))] <- gsub("\"+", "", dat2[, grep("description", tolower(colnames(dat2)))])
	dat2[, grep("description", tolower(colnames(dat2)))] <- gsub("\\#+", "", dat2[, grep("description", tolower(colnames(dat2)))])
}

if(length(grep("symbol", tolower(colnames(dat2)))) > 0) {
	dat2[, grep("symbol", tolower(colnames(dat2)))] <- gsub("\'+", "", dat2[, grep("symbol", tolower(colnames(dat2)))])
	dat2[, grep("symbol", tolower(colnames(dat2)))] <- gsub("\"+", "", dat2[, grep("symbol", tolower(colnames(dat2)))])
	dat2[, grep("symbol", tolower(colnames(dat2)))] <- gsub("\\#+", "", dat2[, grep("symbol", tolower(colnames(dat2)))])
}

# Write data out
write.table(dat2, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)



