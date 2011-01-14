# ANALYSIS Annotation/"Add annotations to data" (Annotates the genes, and adds the results to the datafile. 
# Currently, this function only works with Agilent, Affymetrix and Illumina data.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT data-with-annotations.tsv


# Adds the annotation to the data
# JTT 21.1.2009
#
# MG 25.10.2010
# modified to cope with annomalies in Description names

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
dat<-read.table(file, header=T, sep="\t", row.names=1, quote="", comment.char="")

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

# Fixes an issue with the ' sign appearing in the Description column
# that causes troubles for downstream tools
datannot$Description <- gsub("'", "", datannot$Description)

# Writing out the annotated data
write.table(datannot, file="data-with-annotations.tsv", sep="\t", row.names=T, col.names=T, quote=F)
