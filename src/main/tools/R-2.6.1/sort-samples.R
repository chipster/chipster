# ANALYSIS Utilities/"Sort samples" (Sorts samples according to a phenodata column. The column should contain
# numerical values, since the samples are sorted in ascending order according to the values.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT sort-samples.tsv
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column specifying how to sort)


# Sort samples
# JTT 6.2.2008

# Loads libraries
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Loads data (which file to search)
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Extracts the column from phenodata
groups<-phenodata[,grep(column, colnames(phenodata))]

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Sorting both data and phenodata
dat2<-dat2[,order(groups)]
if(ncol(calls)>0) {
   calls<-calls[,order(groups)] 
}
phenodata2<-phenodata[order(groups),]

# Writing data and phenodata to disk
if(ncol(calls)>0) {
   write.table(data.frame(dat2, calls), file="sort-samples.tsv", sep="\t", row.names=T, col.names=T, quote=F)
} else {
   write.table(data.frame(dat2), file="sort-samples.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}
write.table(phenodata2, file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)
