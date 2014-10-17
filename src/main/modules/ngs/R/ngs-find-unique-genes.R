# TOOL ngs-find-unique-genes.R: "Find unique and annotated genes" (This tool takes a list of ENSEMBL gene identifiers, removes duplicates and fetches annotation information. The output file is compatible with downstream tools used for example for pathway analysis.)
# INPUT ensembl-list.tsv: "Table with ENSEMBL identifiers" TYPE GENERIC 
# OUTPUT unique-genes.tsv: "Table listing the unique genes that can be mapped with gene symbols and entrez gene ids." 
# PARAMETER species: "Genome" TYPE [Human: "Human", Mouse: "Mouse", Rat: "Rat", Zebrafish: "Zebrafish"] DEFAULT Human (The genome to use for fetching annotations.)

# MG, 13.07.2011                                       
# EK, 17.05.2012, fixed parameters
# EK, 22.05.2012, fixed parameters
# Removes duplicate ENSEMBL identifiers, converts them to unique Entrez identifiers, and adds gene annotations to the output 
# MK, 08.05.2012, Added new genomes 
# AMS, 17.10.2014, attribute name change in ensembl: external_gene_id -> external_gene_name

# Set up testing parameters
# species <- "Human"

# Load the libraries
library(ChIPpeakAnno)
library(biomaRt)

# Load the annotation data
if (species == "Human") {
	ensembl_dataset <- "hsapiens_gene_ensembl"
}
if (species == "Mouse") {
	ensembl_dataset <- "mmusculus_gene_ensembl"
}
if (species == "Rat") {
	ensembl_dataset <- "rnorvegicus_gene_ensembl"
}
if (species == "Zebrafish") {
	ensembl_dataset <- "drerio_gene_ensembl"
}

# Read in data and extract the identifiers
#ensembl_list <- read.table (file="ensembl-list.tsv", sep="\t", header=T)
#ensembl_list <- read.table (file="ensembl-list.tsv", sep="\t", header=T, row.names=1, quote="")

ensembl_list <- read.table (file="ensembl-list.tsv", sep="\t", header=T, row.names=NULL, quote="")
rownames(ensembl_list) <- make.names(ensembl_list[,1], unique=T)
ensembl_list <- ensembl_list[,-1]

if("ensembl_id" %in% colnames(ensembl_list)) {
	ensembl_id_list <- as.character(unique(ensembl_list$ensembl_id))
} else {
	ensembl_id_list <- rownames(ensembl_list)
}

# Fetch additional gene info from BioMart database
ensembl_annotations <- useMart("ensembl", dataset=as.character(ensembl_dataset))
#gene_annotations <- getBM(filters="ensembl_gene_id", values=ensembl_id_list, attributes=c("ensembl_gene_id","external_gene_id","description","entrezgene"), mart=ensembl_annotations)
gene_annotations <- getBM(filters="ensembl_gene_id", values=ensembl_id_list, attributes=c("ensembl_gene_id","external_gene_name","description","entrezgene"), mart=ensembl_annotations)

# Keep only the ones that have a unique ensembl gene id and actually map to unique entrez id:s
gene_annotations_na <- na.omit(gene_annotations)
unique_ensembl_id <- unique(levels(as.factor(gene_annotations_na[,1])))
matches <- na.omit(match (unique_ensembl_id, gene_annotations_na[,1], nomatch=NA))
gene_annotations_unique <- gene_annotations_na[matches,]
unique_entrez_id <- unique(levels(as.factor(gene_annotations_unique[,4])))
matches <- na.omit(match (unique_entrez_id, gene_annotations_unique[,4], nomatch=NA))
gene_annotations_unique <- gene_annotations_unique[matches,]
rownames(gene_annotations_unique) <- gene_annotations_unique[,1]
gene_annotations_unique <- gene_annotations_unique[,2:4]
names(gene_annotations_unique)[1] <- "symbol"

# Write output
write.table(gene_annotations_unique, file="unique-genes.tsv", sep="\t", col.names=T, row.names=T, quote=F)

# EOF

