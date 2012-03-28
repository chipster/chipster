# TOOL add-locations-to-data.R: "Add genomic location information to data" (Annotates the genes with information about chromosomal location, and adds the results to the datafile.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT data-with-locations.tsv: data-with-locations.tsv 
# PARAMETER annotate_with: annotate using TYPE [probe_id: probe ID, gene_symbol: gene symbols] DEFAULT gene_symbols (Should the probe identifiers be used to fetch the location of the corresponding gene targets or should gene symbols be used directly, if available.)

# MG 16.3.2012

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
file <- c("normalized.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)

# Extract probes in query list
probes_query <- rownames(dat)

# Load the chromosome number info into memory
lib2 <- sub('.db','',lib)
lib3 <- paste(lib2, "CHR", sep="")
env <- get(lib3)
chr_name <- unlist(mget(probes_query, env))

# Load the chromosome start position info into memory
lib3 <- paste(lib2, "CHRLOC", sep="")
env <- get(lib3)
chr_start <- unlist(mget(probes_query, env))

# Load the chromosome end position info into memory
lib3 <- paste(lib2, "CHRLOCEND", sep="")
env <- get(lib3)
chr_end <- unlist(mget(probes_query, env))

# Figure out the strand information and remove sign from coordinates
is_negative <- ifelse ((as.numeric(chr_start)-1)>-1, "-", "+")
chr_start <- gsub("-", "", x=as.character(chr_start))
chr_end <- gsub("-", "", x=as.character(chr_end))
length_info <- as.numeric (chr_end)-as.numeric(chr_start)+1

# Merge things together and output file
output_file <- data.frame(dat, chromosome=chr_name, start=chr_start, end=chr_end, length=length_info, strand=is_negative)
write.table(datannot, file="data-with-annotations.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# EOF
