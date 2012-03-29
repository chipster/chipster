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
library(biomaRt)

# Loads the data
file <- c("normalized.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)

# Extract probes in query list
probes_query <- rownames(dat)

# Load the gen symbol info into memory
lib2 <- sub('.db','',lib)
lib3 <- paste(lib2, "SYMBOL", sep="")
env <- get(lib3)
gene_symbols <- unlist(mget(probes_query, env))


# Fetch the gene symbols and descriptions
if (species=="human") {
	dataset <- "hsapiens_gene_ensembl"
}
if (species=="mouse") {
	dataset <- "mmusculus_gene_ensembl"
}
if (species=="rat") {
	dataset <- "rnorvegicus_gene_ensembl"
}
ensembl <- useMart("ensembl", dataset=dataset)
annotated_genes <- getBM(mart=ensembl, attributes=c("hgnc_symbol","chromosome_name","start_position","end_position"), filters="hgnc_symbol", values=gene_symbols)
# remove chromosome entries containing the "_" character
annotated_genes <- annotated_genes [-grep( pattern="_", annotated_genes$chromosome_name),]
rownames(annotated_genes) <- annotated_genes$hgnc_symbol

# Match the list of input gene ids with the annotations
chr_name <- character(length(probes_query))
chr_start <- character(length(probes_query))
chr_end <- character(length(probes_query))
result_table <- data.frame(chromosome=chr_name, start=chr_start, end=chr_end, dat, stringsAsFactors = FALSE)

for (gene_count in 1:length(probes_query)) {
	chr_name[gene_count] <- annotated_genes[gene_symbols[gene_count],2]
	chr_start[gene_count] <- annotated_genes[gene_symbols[gene_count],3]
	chr_end[gene_count] <- annotated_genes[gene_symbols[gene_count],4]
#	result_table [gene_id==gene_entrezid,3] <- gene_symbol
#	result_table [gene_id==gene_entrezid,4] <- gene_description
}	
result_table <- data.frame(chromosome=chr_name, start=chr_start, end=chr_end, dat, stringsAsFactors = FALSE)

names(result_table) <- c("miRNA","entrez_id","symbol","description")

# add rownames to allow use of Venn diagrams
concatenate.if.not.equal <- function(x) {
	x <- unique(x)
	paste(x, collapse=';')
}
result_table$pair <- paste(result_table$miRNA, result_table$entrez_id, sep=';')
result_table2 <- aggregate(result_table[,-5], list(result_table[,5]), concatenate.if.not.equal)
result_table2$count <- aggregate(result_table[,5], list(result_table[,5]), length)$x
rownames(result_table2) <- result_table2[,1]
result_table2[,1] <- NULL


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
names(chr_start) <- names(chr_name)
names(chr_end) <- names(chr_name)
names(length_info) <- names(chr_name)

# Collapse multiple gene locations into one
probe_names <- rownames(dat)
number_probes <- length(probe_names)
chr_name_collapsed <- numeric(number_probes)
for (probe_count in 1:number_probes) {
	chr_start_subset <- min(as.numeric(chr_start[probe_names[probe_count]])) 
	chr_name_collapsed[probe_count] <- chr_start_subset
}

# Merge things together and output file
output_file <- data.frame(dat, chromosome=chr_name, start=chr_start, end=chr_end, length=length_info, strand=is_negative)
write.table(datannot, file="data-with-annotations.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# EOF
