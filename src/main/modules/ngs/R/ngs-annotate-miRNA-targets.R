# TOOL ngs-annotate-miRNA-targets.R: "Find miRNA targets" (Fetches the predicted gene targets of a list of miRNA ID:s in one of the following databases: miRanda, miRbase, miRtarget2, PicTar, TarBase or targetScan.)
# INPUT normalized_mirna.tsv: normalized_mirna.tsv TYPE GENERIC 
# OUTPUT mirna_targets.tsv: mirna_targets.tsv 
# PARAMETER species: Species TYPE [human: human, mouse: mouse, rat: rat] DEFAULT human (The species for which the miRNA:s have been analyzed.)
# PARAMETER database: Database TYPE [miranda: miRanda, mirbase: miRbase, mirtarget2: miRtarget2, pictar: PicTar, tarbase: TarBase, targetscan: targetScan] DEFAULT mirtarget2 (The database from which to search for predicted target genes.)

# MG, 22.2.2010
# IS, 13.10.2010, added rownames of miRNA-target-pairs to allow intersecting through Venn diagrams
# MK, 26.04.2013, supports now files having IDs in ID-field and/or having alignment information in their ID-field
# MK, 12.04.2013, Added new organisms

# Load the required libraries
library(biomaRt)

# Reads the data
mirna_data <- read.table("normalized_mirna.tsv", sep="\t", header=T)

# Extracts the identifiers
if("id" %in% colnames(mirna_data)) {
	mirna_id <- as.character(mirna_data$id)
} else {
	mirna_id <- as.character(rownames(mirna_data))
}

# If convert genomic BAM file has been used, table has a column which name is sequence
if("sequence" %in% colnames(mirna_data) && (length(grep(mirna_data$sequence[1], mirna_id[1])) > 0)) {
	mirna_id_list <- strsplit(as.vector(mirna_id), "_")
	mirna_id <- NULL
	for(i in 1:length(mirna_id_list)) {
		#remove last three sections of the id
		mirna_id <- c(mirna_id, paste(unlist(mirna_id_list[i])[1:((length(unlist(mirna_id_list[i])))-3)], collapse="_"))
	}
}

# Fetch miRNA targets
if(species == "human") {
	dataset <- "hsapiens_gene_ensembl"
	library(RmiR.Hs.miRNA)
	dbconn <- RmiR.Hs.miRNA_dbconn()
	full_database <- dbReadTable(dbconn, database)[, 1:2]
	target_genes <- full_database[full_database$mature_miRNA %in% mirna_id,]

} else if (species=="mouse") {
	if (database != "targetscan") {
		stop('CHIPSTER-NOTE: Only TargetScan database is available for mouse')
	}

	dataset <- "mmusculus_gene_ensembl"
	library(targetscan.Hs.eg.db)

	target_genes <- NULL
	for(i in 1:length(mirna_id)) {
		mirna.fam <- try(mget(as.vector(mirna_id[i]), revmap(targetscan.Hs.egFAMILY2MIRBASE)), silent=T)
		if(class(mirna.fam) != "try-error") {
			mirna.targets <- try(mget(unlist(mirna.fam), revmap(targetscan.Hs.egTARGETS)), silent=T)
			if(class(mirna.targets) != "try-error") {	
				mirna.targets <- unique(unlist(mirna.targets))
				if(!(is.na(mirna.targets))) {
					target_genes <- rbind(target_genes, as.data.frame(cbind(mature_miRNA=rep(mirna_id[i], length(mirna.targets)), gene_id=mirna.targets)))
				}
			}
		}
	}
} else if (species=="rat") {
	dataset <- "rnorvegicus_gene_ensembl"

	#no genuine rat databases
	library(org.Rn.eg.db)
	ensembl.to.entrez <- as.list(org.Rn.egENSEMBLTRANS2EG)

	targets <- read.table(file.path(chipster.tools.path, "miRNA_mappings", "mirna_mappings_rnorvegicus.txt"), sep="\t")
	targets2 <- as.data.frame(cbind(gene_id=as.numeric(as.character(ensembl.to.entrez [c(as.vector(targets$tran))])), mature_miRNA=as.character(targets[,2])))
	targets2 <- targets2[which(!is.na(targets2$gene_id)),]
	target_genes <- targets2[tolower(targets2$mature_miRNA) %in% mirna_id, ]
}

gene_id <- target_genes$gene_id
mirna_id_2 <- target_genes$mature_miRNA

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
#most likely, number of predicted binding sites
result_table2$sites <- aggregate(result_table[,5], list(result_table[,5]), length)$x
rownames(result_table2) <- result_table2[,1]
result_table2[,1] <- NULL

if (species=="mouse") {
	result_table2$sites <- rep(NA, nrow(result_table2))
}
if (species=="rat") {
	result_table2$sites <- rep(NA, nrow(result_table2))
}

# Write out the output file
write.table(result_table2, file="mirna_targets.tsv", sep="\t", quote=FALSE)

# EOF
