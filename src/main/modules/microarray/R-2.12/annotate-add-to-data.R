# TOOL annotate-add-to-data.R: "Add annotations to data" (Annotates the genes, and adds the results to the datafile. Currently, this function only works with Agilent, Affymetrix and Illumina data.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT data-with-annotations.tsv: data-with-annotations.tsv 
# PARAMETER conditional: conditional TYPE [yes: yes, no: no] DEFAULT no (When run in Conditional mode, each gene is not only annotated to the most specific GO term but also to all its parent terms up to the top of the hierarchy. This option is recommended in conjunction to the tool Hypergeometric test for GO run in the Conditional testing mode.)

# Adds the annotation to the data
# JTT 21.1.2009
#
# MG 25.10.2010
# modified to cope with annomalies in Description names
# MG 25.1.2012
# modified to allow for conditional mode for GO terms

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

# If Conditional mode run this
if (conditional == "yes") {

	# Extract probes in query list
	probes.query <- rownames(dat)

	# Extract the mapping info
	lib2<-sub('.db','',lib)
	lib3 <- paste(lib2, "GO2ALLPROBES", sep="")
	env <- get(lib3)
	probes.2.go <- revmap(as.list(env))[probes.query] 

	# Input the GO ids that correspond to each probe id
	for (count in 1:length(probes.query)) {
		# merge GO id and term definition
		go.list <- unique(probes.2.go[[count]])
		if (is.null(go.list)) go.list <- ""
		go.list <- paste(go.list, Term(go.list))
		# concatenate GO list into single character vector
		go.list <- paste(go.list, collapse="; ")
		datannot$Gene.Ontology[count] <- go.list
	}

}

# Fixes an issue with the ' sign appearing in the Description column
# that causes troubles for downstream tools
datannot$Description <- gsub("'", "", datannot$Description)

# Remove any NA entries and replace with empty string in case there is no match
# between a probe id and GO id
datannot$Gene.Ontology[datannot$Gene.Ontology == " NA"] <- ""

# Writing out the annotated data
write.table(datannot, file="data-with-annotations.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# EOF
