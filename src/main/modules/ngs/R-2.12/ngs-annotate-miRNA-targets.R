# TOOL ngs-annotate-miRNA-targets.R: "Find miRNA targets" (Fetches the predicted gene targets of a list of miRNA ID:s in one of the following databases: miRanda, miRbase, miRtarget2, PicTar, TarBase or targetScan.)
# INPUT normalized_mirna.tsv: normalized_mirna.tsv TYPE GENE_EXPRS 
# OUTPUT mirna_targets.tsv: mirna_targets.tsv 
# PARAMETER species: species TYPE [human: human] DEFAULT human (The species for which the miRNA:s have been analyzed. Mouse and rat will be supported as soon as the RmiR package from Bioconductor supports this.)
# PARAMETER database: database TYPE [miranda: miranda, mirbase: mirbase, mirtarget2: mirtarget2, pictar: pictar, tarbase: tarbase, targetscan: targetscan] DEFAULT mirtarget2 (The database from which to search for predicted target genes.)

# miRNA hypergeometric test for GO
# MG, 22.2.2010
# IS, 13.10.2010, added rownames of miRNA-target-pairs to allow intersecting through Venn diagrams

# Load the required libraries
library(RmiR.Hs.miRNA)
library(biomaRt)

# Reads the data
mirna_data <- read.table("normalized_mirna.tsv", sep="\t", header=T)

# Extracts the identifiers
mirna_id <- as.character(rownames(mirna_data))

# Fetch the complete set of target genes for a specified database
full_database <- dbReadTable(RmiR.Hs.miRNA_dbconn(), database)[, 1:2]
target_genes <- full_database[full_database$mature_miRNA %in% mirna_id,]
gene_id <- target_genes$gene_id
mirna_id_2 <- target_genes$mature_miRNA

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
annotated_genes <- getBM(mart=ensembl, attributes=c("entrezgene","hgnc_symbol","description"), filters="entrezgene", values=gene_id)

# Match the list of transcripts with the annotations
gene_entrezid <- character(length(gene_id))
gene_symbol <- character(length(gene_id))
gene_description <- character(length(gene_id))
result_table <- data.frame(mirna_id_2, gene_id, gene_symbol, gene_description, stringsAsFactors = FALSE)
for (gene_count in 1:length(annotated_genes [,1])) {
	gene_entrezid <- annotated_genes[gene_count,1]
	gene_symbol <- annotated_genes[gene_count,2]
	gene_description <- annotated_genes[gene_count,3] 
	result_table [gene_id==gene_entrezid,3] <- gene_symbol
	result_table [gene_id==gene_entrezid,4] <- gene_description
}	
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

# Write out the output file
write.table(result_table2, file="mirna_targets.tsv", sep="\t", quote=FALSE)

# EOF
