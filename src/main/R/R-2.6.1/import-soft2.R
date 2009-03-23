# ANALYSIS Utilities/"Import from GEO" (Imports a SOFT-formatted datafile directly from GEO. 
# Be sure to specify the chiptype as an Affymetrix chip name, or either Illumina or cDNA.)
# OUTPUT normalized.tsv, phenodata.tsv 
# PARAMETER GDS.name STRING DEFAULT GDS100 (GDS or GSE number of the experiment.)
# PARAMETER chiptype STRING DEFAULT cDNA (Chiptype)
# PARAMETER log2.transform [yes, no] DEFAULT yes (log2 transform the data)


# JTT 9.8.2007

# Loads the libraries
library(GEOquery)

# Loads and parses the SOFT file
if(GDS.name!="empty") {
   gds<-getGEO(GDS.name)
   if(log2.transform=="yes") {
      eset<-GDS2eSet(gds, do.log2=T)
   } else {
      eset<-GDS2eSet(gds, do.log2=F)
   }      
   dat<-exprs(eset)
   colnames(dat)<-paste("chip.", colnames(dat), sep="")
}

if(GDS.name=="empty") {
   stop("You need to specify a valid GDS!")
}

# Writes out a phenodata
sample<-colnames(exprs(eset))
group<-c(rep("", nrow(pData(eset))))
training<-c(rep("", nrow(pData(eset))))
time<-c(rep("", nrow(pData(eset))))
random<-c(rep("", nrow(pData(eset))))
write.table(data.frame(sample=sample, group=group, training=training, chiptype=chiptype), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Writes out a normalized datafile
write.table(as.data.frame(dat), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
