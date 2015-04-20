# TOOL annotate-miRNA-targets.R: "Find miRNA targets" (Fetches the predicted gene targets of a list of miRNA ID:s in one of the following databases: miRanda, miRbase, miRtarget2, PicTar, TarBase or targetScan.)
# INPUT normalized_mirna.tsv: normalized_mirna.tsv TYPE GENE_EXPRS 
# OUTPUT mirna_targets.tsv: mirna_targets.tsv 
# PARAMETER species: Species TYPE [human: human] DEFAULT human (The species for which the miRNA:s have been analyzed. Mouse and rat will be supported as soon as the RmiR package from Bioconductor supports this.)
# PARAMETER database: Database TYPE [miranda: miRanda, mirbase: miRBase, mirtarget2: miRtarget2, pictar: PicTar, tarbase: TarBase, targetscan: targetScan] DEFAULT mirtarget2 (The database from which to search for predicted target genes.)

# MG, 22.02.2010
# IS, 13.10.2010, added rownames of miRNA-target pairs to allow intersecting through Venn diagrams
# MK, 09.12.2013, performance of the script improved. Because of this, information associated with the first matching row is used instead of the last matching row 

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

# Match the list of genes with the annotations
result_table <- data.frame(cbind(mirna_id_2, gene_id, matrix("", nrow=length(gene_id), ncol=2)), stringsAsFactors = FALSE)
names(result_table) <- c("miRNA","entrez_id","symbol","description")
# match function reports the first match in annotated genes, not the last as done previously 
result_table[, c(3,4)] <- annotated_genes[match(gene_id, annotated_genes[,1]), 2:3]

if(length(which(is.na(result_table[,3]))) > 0) {
	result_table[which(is.na(result_table[,3])),3] <- ""
	result_table[which(is.na(result_table[,4])),4] <- ""
}

# add rownames to allow use of Venn diagrams
concatenate.if.not.equal <- function(x) {
  x <- unique(x)
  paste(x, collapse=';')
}

result_table$pair <- paste(result_table$miRNA, result_table$entrez_id, sep=';')
result_table2 <- aggregate(result_table[,-5], list(result_table[,5]), concatenate.if.not.equal)

result_counts <- table(result_table[,5]);
result_table2$count <- rep(0, nrow(result_table2))
result_table2$count <- as.vector(result_counts[match(result_table2[,1], names(result_counts))])

#result_table2$count <- aggregate(result_table[,5], list(result_table[,5]), length)$x
rownames(result_table2) <- result_table2[,1]
result_table2[,1] <- NULL

# Write out the output file
write.table(result_table2, file="mirna_targets.tsv", sep="\t", quote=FALSE)

# EOF
