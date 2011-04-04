# ANALYSIS Utilities/"Extract genes from KEGG pathway" (Fetches the genes that belong to a given KEGG pathway, defined either by ID or description.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv
# OUTPUT extracted-from-KEGG.tsv
# PARAMETER match.term STRING DEFAULT empty (The identifier for the KEGG pathway to extract genes from.)

# Extract the rows from a table that contain probes that map genes annotated
# to a given KEGG pathway
#
# MG 26.11.2010


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
lib2<-sub('.db','',lib)
library(package=lib, character.only=T)

# Loads the data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Extract the mapping info
lib3 <- paste(lib2, "PATH2PROBE", sep="")
env <- get(lib3)
kegg.2.probes <- as.list(env) 

# Extract all probes belonging to query
probes.list <- kegg.2.probes[match.term]
probes.kegg <- unlist(probes.list)

# Extract probes in query list
probes.query <- rownames(dat)

# Extract common probes
match.indices <- match (probes.kegg, probes.query, nomatch=0)
match.indices <- unique (match.indices[match.indices>0])

# Extract data
dat2 <- dat[match.indices,]

# Writing out the extracted data
write.table(dat2, file="extracted-from-KEGG.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# EOF