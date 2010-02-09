# ANALYSIS Utilities/"Import from ArrayExpress" (Import Affymetrix raw data directly from ArrayExpress, and normalise it using the RMA algorithm.)
# OUTPUT normalized.tsv, phenodata.tsv 
# PARAMETER accession STRING DEFAULT E-MEXP-1422 (Accession number of the experiment.)


# 13.1.2010
# JTT

# Loads the libraries
library(ArrayExpress)

# Loads the data
dat<-ArrayExpress(accession)

# Normalizes the raw data
if(class(dat)=="AffyBatch") {
   library(affy)
   dat2<-rma(dat)
   dat2<-exprs(dat2)
   colnames(dat2)<-paste("chip.", colnames(dat2), sep="")
   chiptype<-dat@annotation
}

# Writing out data
a<-try(library(paste(chiptype, ".db", sep=""), character.only=T))
if(chiptype!="empty" & class(a)!="try-error") {
   # Including gene names to data
   lib2<-sub('.db','',chiptype)
   symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep="")))))[rownames(dat2),])
   genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep="")))))[rownames(dat2),])
   # Writes the results into a file
   write.table(data.frame(symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
} 

if(chiptype=="empty" | class(a)=="try-error") {
   write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
}

# Writes out a phenodata table
sample<-rownames(pData(dat))
group<-c(rep("", nrow(pData(dat))))
chiptype<-paste(chiptype, ".db", sep="")
write.table(data.frame(sample=sample, chiptype=chiptype, group=group, pData(dat)), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

