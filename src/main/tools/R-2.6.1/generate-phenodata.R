# ANALYSIS Utilities/"Generate phenodata" (If run on a prenormalized file, generates a blank phenodata for it.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT phenodata.tsv, normalized.tsv
# PARAMETER chiptype STRING DEFAULT empty (Chiptype)


# Combines two different tables using gene names
# JTT 6.7.2006

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-data.frame(dat[,grep("chip", names(dat))])

# Generates the variables
sample<-colnames(dat2)
group<-c(rep("", ncol(dat2)))
if(chiptype=="empty") {
   chiptype<-c("empty")
} else {
   chiptype<-chiptype
}

# Writes out the data and the phenodata table
write.table(dat, file="normalized.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)
