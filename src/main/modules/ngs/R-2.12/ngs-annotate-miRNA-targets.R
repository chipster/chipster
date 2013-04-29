# TOOL ngs-annotate-miRNA-targets.R: "Find miRNA targets" (Fetches the predicted gene targets of a list of miRNA ID:s in one of the following databases: miRanda, miRbase, miRtarget2, PicTar, TarBase or targetScan.)
# INPUT normalized_mirna.tsv: normalized_mirna.tsv TYPE GENERIC 
# OUTPUT mirna_targets.tsv: mirna_targets.tsv 
# PARAMETER species: species TYPE [human: human] DEFAULT human (The species for which the miRNA:s have been analyzed. Mouse and rat will be supported as soon as the RmiR package from Bioconductor supports this.)
# PARAMETER database: database TYPE [miranda: miranda, mirbase: mirbase, mirtarget2: mirtarget2, pictar: pictar, tarbase: tarbase, targetscan: targetscan] DEFAULT mirtarget2 (The database from which to search for predicted target genes.)

# MG, 22.2.2010
# IS, 13.10.2010, added rownames of miRNA-target-pairs to allow intersecting through Venn diagrams
# MK, 26.04.2013, supports now files having IDs in ID-field and/or having alignment information in their ID-field

# Load the required libraries
library(RmiR.Hs.miRNA)
library(biomaRt)

# Reads the data
mirna_data <- read.table("normalized_mirna.tsv", sep="\t", header=T)

# Extracts the identifiers

#if has ID, field try that
if("id" %in% colnames(mirna_data)) {
	mirna_id <- mirna_data$id
} else {
	mirna_id <- as.character(rownames(mirna_data))
}
	
# Fetch the complete set of target genes for a specified database
full_database <- dbReadTable(RmiR.Hs.miRNA_dbconn(), database)[, 1:2]
target_genes <- full_database[full_database$mature_miRNA %in% mirna_id,]

# If convert genomic BAM file has been used, IDs contain read information preveting direct matching. To correct this,
# database IDs are first compared to user IDs (allowing partial matches), after which those user IDs having more than
# one match are rejected as I cannot be known which of the many is the correct
if(nrow(target_genes) == 0) {
	mature_miRNA_mod 		<- unique(paste(full_database$mature_miRNA, "_", sep=""))
	mirna_id_mat 	  		<- data.frame(counts=rep(0, length(mirna_id)), names=rep("", length(mirna_id)), stringsAsFactors = FALSE)
	rownames(mirna_id_mat) 	<- mirna_id

	for(i in 1:length(mature_miRNA_mod)) {
		grep.positions <- grep(paste("^",mature_miRNA_mod[i], sep=""), mirna_id)
		if(length(grep.positions) > 0) {
			mirna_id_mat[grep.positions, 1] <- mirna_id_mat[grep.positions, 1] + 1
			mirna_id_mat[grep.positions, 2] <- paste(mirna_id_mat[grep.positions, 2], mature_miRNA_mod[i], sep="")		
		}
	}
	target_genes <- full_database[paste(full_database$mature_miRNA, "_", sep="") %in% mirna_id_mat[mirna_id_mat[,1] >= 1, 2],]
}
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

if(length(gene_id) < 1) {
	stop("CHIPSTER-NOTE: Your analysis was terminated, because no matches were found between your and RmiR database miRNA identifiers")
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
