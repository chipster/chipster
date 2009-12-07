# ANALYSIS Pathways/"KEGG enrichment for miRNA targets" (Performs a statistical test for enrichmens of
# KEGG pathways in the predicted gene targets of a list of miRNA ID:s.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT hyperg_kegg.tsv
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.01 (P-value threshold)
# PARAMETER species [human] DEFAULT human (the species for which the miRNA:s have been analyzed)

# POSSIBLE summary.feature [gene, transcript] DEFAULT gene (should the targets for the miRNA:s be transcripts or genes?)


# miRNA hypergeometric test for KEGG
# MG, 4.11.2009

# force "transcript" mode
summary.feature <- "transcript"

# Reads the data
dat<-read.table("normalized.tsv", sep="\t", header=T)

# Extracts the identifiers
id<-as.character(rownames(dat))

# Translate parameter settings
if (species=="human" & summary.feature=="gene") {
	dataset <- "hsapiens_gene_ensembl"
}
if (species=="human" & summary.feature=="transcript") {
	dataset <- "hsapiens_transcript_ensembl"
}
if (species=="mouse" & summary.feature=="gene") {
	dataset <- "mmusculus_gene_ensembl"
}
if (species=="mouse" & summary.feature=="transcript") {
	dataset <- "mmusculus_transcript_ensembl"
}
if (species=="rat" & summary.feature=="gene") {
	dataset <- "rnorvegicus_gene_ensembl"
}
if (species=="rat" & summary.feature=="transcript") {
	dataset <- "rnorvegicus_transcript_ensembl"
}

# Read in the CORNA library, which contains the functions to map miRNA:s to targets
# and performs hypergeometric test for enrichment of GO terms
library(CORNA)

# Download the mapping of miRNA to its targets from Sanger institute
#if (species=="human") {
#	targets <- miRBase2df.fun(url="ftp://ftp.sanger.ac.uk/pub/mirbase/targets/v5/arch.v5.txt.homo_sapiens.zip")
#}
#if (species=="mouse") {
#	targets <- miRBase2df.fun(url="ftp://ftp.sanger.ac.uk/pub/mirbase/targets/v5/arch.v5.txt.mus_musculus.zip")
#}
#if (species=="rat") {
#	targets <- miRBase2df.fun(url="ftp://ftp.sanger.ac.uk/pub/mirbase/targets/v5/arch.v5.txt.rattus_norvegious.zip")
#}


# Download the mapping of miRNA to its targets from locally installed files
path.mappings <- c(file.path(chipster.tools.path, "miRNA_mappings"))
if (species=="human") {
	targets <- read.table(file.path(path.mappings, "mirna_mappings_hsapiens.txt"), sep="\t")
}
if (species=="mouse") {
	targets <- read.table(file.path(path.mappings, "mirna_mappings_mmusculus.txt"), sep="\t")
}
if (species=="rat") {
	targets <- read.table(file.path(path.mappings, "mirna_mappings_rnorvegious.txt"), sep="\t")
}


# obtain a link from transcript to gene from biomart
#tran2gene <- BioMart2df.fun(
#		biomart="ensembl",
#		dataset=dataset,
#		col.old=c("ensembl_transcript_id",
#				"ensembl_gene_id"),
#		col.new=c("tran", "gene"))

# obtain a link from transcript to gene from locally installed files
#tran2gene <- read.table(file="transcripts_to_genes_hsapiens.txt", sep="\t")
tran2gene <- read.table(file=file.path(path.mappings,"transcripts_to_genes_hsapiens.txt"), sep="\t")

# link microRNAs to genes instead of transcripts
mir2gene <- corna.map.fun(targets,
		tran2gene,
		"gene",
		"mir")

# get those genes associated with a list of regulated miRNA:s
sample.list <- corna.map.fun(mir2gene,
		id,
		"mir",
		"gene")

# read pathway information from KEGG for homo sapiens
#gene2path <- KEGG2df.fun(org="hsa")

# read pathway information from local file for homo sapiens
#gene2path <- read.table(file="genes_to_kegg_hsapiens.txt", sep="\t")
gene2path <- read.table(file=file.path(path.mappings,"genes_to_kegg_hsapiens.txt"), sep="\t")

# define the sample of genes to test, i.e. only those gene targets that have pathway data
sample.list <- intersect(sample.list, unique(gene2path$gene))

# create data frame of pathway ids and names
path2name <- unique(gene2path[c("path", "name")])

# perform the test
test <- corna.test.fun(
		sample.list,
		unique(gene2path$gene),
		gene2path,
		hypergeometric=T,
		fisher=T,
		chi.square=T,
		#fisher.alternative="greater",
		min.pop=10,
		sort="fisher",
		desc=path2name)

# fetch significant pathways
significant.kegg <- test[test$hypergeometric<=0.05,]

# write results table
write.table(significant.kegg, file="hyperg_kegg.tsv", sep="\t", quote=F)


