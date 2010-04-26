# ANALYSIS Utilities/"Extract samples from dataset" (Extracts samples from a dataset. Saves the extracted samples as a
# new dataset. The samples to be extracted are coded with 1 in one column of the phenodata. The samples to be deleted
# are coded with 0 in the same column.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT extract.tsv, phenodata.tsv
# PARAMETER column.extract METACOLUMN_SEL DEFAULT group (Phenodata column containing the samples to be extracted)

# Extracts subset of samples from a dataset
# JTT 19.10.2007
#
# modified, MG, 20.4.2010 to include annotation info

# Loads the data file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Check if there is annotation info available and if so extract it
annotations <- dat[,-c(grep("chip",names(dat)), grep("flag", names(dat)))]
if (length(annotations)>0) {
	rownames(annotations) <- rownames(dat)
}

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Extract the data from the phenodata column
extract<-phenodata[,pmatch(column.extract,colnames(phenodata))]

# extract<-phenodata[,grep(column.extract, colnames(phenodata))]

# Sanity checks
if(length(unique(extract))>2) {
   stop("CHIPSTER-NOTE: You have specified more than two groups! You need to define exactly two groups.")
}
if(max(extract>1)) {
   stop("CHIPSTER-NOTE: The groups should be defined with values of 0 and 1! You have numbers larger than 1 in the definitions.")
}

# Extracting the samples
dat2<-dat2[,which(extract==1)]
if(ncol(calls)>=1) {
   calls2<-calls[,which(extract==1)]
}
phenodata2<-phenodata[which(extract==1),]

# Writing the data to disk
if (length(annotations)>0) {
	dat2 <- data.frame(annotations,dat2)
}
if(ncol(calls)>=1) {
   write.table(data.frame(dat2, calls2), file="extract.tsv", sep="\t", row.names=T, col.names=T, quote=F)
} else {
   write.table(data.frame(dat2), file="extract.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}
write.table(phenodata2, file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)
