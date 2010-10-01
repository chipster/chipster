# ANALYSIS Utilities/"Extract genes from GO term" (Fetches the genes that belong to a given GO term.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT data-with-annotations.tsv
# PARAMETER match.term STRING DEFAULT empty (String to search for.)


# Fetches genes for a given GO term
# MG 9.8.2010

# Reads the chiptype from phenodata table
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
if(phenodata$chiptype[1]!="cDNA" | phenodata$chiptype[1]!="Illumina") {
   # Saves the chiptype into object lib
   lib<-phenodata$chiptype[1]
   lib<-as.character(lib)
}

# Account for the fact that annotation packages are from version 2.3 of Bioconductor
# named with an ".db" suffix. Add the suffix when missing to support data files
# from Chipster 1.3 and earlier. 
if (length(grep(".db", lib)) == 0 & length(grep("pmcdf", lib)) == 0) {
  lib <- paste(lib, ".db", sep="")
}

# Loads the correct annotation library
library(package=lib, character.only=T)
library(annaffy)

# Loads the data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Extract the mapping info
translation_map <- 

# Creating annotations from the library
annot.cols<-aaf.handler()
annot.table<-aafTableAnn(row.names(dat), lib, annot.cols)
saveText(annot.table, "annotations.tsv")

# Merging annotations and data
annot<-read.table("annotations.tsv", header=T, sep="\t", row.names=1, quote="")
annot$Pathway<-gsub("\'", "", annot$Pathway)
annot$Gene.Ontology<-gsub("\'", "", annot$Gene.Ontology)
datannot<-merge(dat, annot, by.x="row.names", by.y="row.names")
rownames(datannot)<-datannot[,1]
datannot<-datannot[,-1]

# Writing out the annotated data
write.table(datannot, file="data-with-annotations.tsv", sep="\t", row.names=T, col.names=T, quote=F)
