# ANALYSIS Utilities/"Sort samples" (Sorts samples according to a phenodata column. The column should contain
# numerical values, since the samples are sorted in ascending order according to the values.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT sort-samples.tsv, phenodata-sorted.tsv
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column specifying how to sort)

# Sort samples
# JTT 6.2.2008
#
# MG, 16.11.2010
# modified to also generate a re-ordered phenodata file to reflect the re-ordered data

# Default parameters
#column<-"group"
 
# Loads libraries
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Loads data (which file to search)
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Extracts the column from phenodata
groups<-phenodata[,pmatch(column,colnames(phenodata))]

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Get indices
indices_calls <- grep("flag", names(dat))
indices_dat <- grep("chip", names(dat))

# Sorting both data and phenodata
dat2<-dat2[,order(groups)]
if(ncol(calls)>0) {
	calls<-calls[,order(groups)] 
}
phenodata2<-phenodata[order(groups),]

# Fill in the ordered data and flag values
dat3 <- dat
dat3[,indices_dat] <- dat2
names(dat3) [indices_dat] <- names(dat2)
if(ncol(calls)>0) {
	dat3[,indices_calls] <- calls
	names(dat3) [indices_calls] <- names(calls)
}

# Writing data and phenodata to disk
#if(ncol(calls)>0) {
#   write.table(data.frame(dat2, calls), file="sort-samples.tsv", sep="\t", row.names=T, col.names=T, quote=F)
#} else {
#   write.table(data.frame(dat2), file="sort-samples.tsv", sep="\t", row.names=T, col.names=T, quote=F)
#}
write.table(dat3, file="sort-samples.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(phenodata2, file="phenodata-sorted.tsv", sep="\t", row.names=F, col.names=T, quote=F)
