# TOOL ngs-find-unique-genes.R: "Find unique and annotated genes" (This tool takes a list of ENSEMBL gene identifiers, removes duplicates and fetches annotation information. The output file is compatible with downstream tools used for example for pathway analysis.)
# INPUT ensembl-list.tsv: "Table with ENSEMBL identifiers" TYPE GENERIC 
# OUTPUT unique-genes.tsv: "Table listing the unique genes that can be mapped with gene symbols and entrez gene ids." 
# PARAMETER species: Species TYPE [Human, Mouse, Rat] DEFAULT Human (The species of the genome to use for fetching annotationsan.)

#######################################################
#                                                     #
# MG, 13.7.2011                                       #
#                                                     #
# Tool that removes duplicate genes from a list of    #
# ENSEMBL identifiers, coverts them to unique entrez  #
# identifiers and adds gene annotations to the output #
#                                                     #
#######################################################

# Set up testing parameters
# species <- "Human"

# Load the libraries
library(ChIPpeakAnno)
library(biomaRt)

# Load the annotation data
if (species == "Human") {
	data(TSS.human.GRCh37)
	annotations <- TSS.human.GRCh37
	ensembl_dataset <- "hsapiens_gene_ensembl"
	filter <- "ens_hs_gene"
}
if (species == "Nouse") {
	data(TSS.mouse.NCBIM37)
	annotations <- TSS.mouse.NCBIM37
	ensembl_dataset <- "mmusculus_gene_ensembl"
	filter <- "ens_mm_gene"
}
if (species == "Rat") {
	data(TSS.rat.RGSC3.4)
	annotations <- TSS.rat.RGSC3.4
	ensembl_dataset <- "rnorvegicus_gene_ensembl"
	filter <- "ens_rn_gene"
}

# Read in data and extract the identifiers
ensembl_list <- read.table (file="ensembl-list.tsv", sep="\t", header=T)
ensembl_id_list <- as.character(unique(ensembl_list$ensembl_id))

# Fetch additional gene info from BioMart database
ensembl_annotations <- useMart("ensembl", dataset=as.character(ensembl_dataset))
gene_annotations <- getBM(filters="ensembl_gene_id", values=ensembl_id_list, attributes=c("ensembl_gene_id","hgnc_symbol","description","entrezgene"), mart=ensembl_annotations)

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

