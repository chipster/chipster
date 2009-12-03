# ANALYSIS Pathways/"miRNA, enriched GO" (Performs a statistical test for enrichments of
# GO terms in the predicted gene targets of a list of miRNA ID:s.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT hyperg_go_mf.tsv, hyperg_go_bp.tsv hyperg_go_cc.tsv
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value threshold)
# PARAMETER species [human, mouse, rat] DEFAULT human (the species for which the miRNA:s have been analyzed)

# POSSIBLE summary.feature [gene, transcript] DEFAULT gene (should the targets for the miRNA:s be transcripts or genes?)


# miRNA hypergeometric test for GO
# MG, 3.10.2009

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

# Fetch the predicted targets for the list of miRNA:s
targets.list <- corna.map.fun(targets, id, "mir", "tran", all=TRUE)

# Get links between targets and GO terms from BioMart:
#if (summary.feature=="gene") {
#	tran2gomf  <- BioMart2df.fun(biomart="ensembl",  
#			dataset=dataset,
#			col.old=c("ensembl_gene_id", "go_molecular_function_id"),
#			col.new=c("gene", "gomf"))
#	tran2gobp  <- BioMart2df.fun(biomart="ensembl",  
#			dataset=dataset,
#			col.old=c("ensembl_gene_id", "go_biological_pathway_id"),
#			col.new=c("gene", "gobp"))
#	tran2gocc  <- BioMart2df.fun(biomart="ensembl",  
#			dataset=dataset,
#			col.old=c("ensembl_gene_id", "go_cellular_component_id"),
#			col.new=c("gene", "gocc"))
#}
#if (summary.feature=="transcript") {
#	tran2gomf  <- BioMart2df.fun(biomart="ensembl",  
#			dataset=dataset,
#			col.old=c("ensembl_transcript_id", "go_molecular_function_id"),
#			col.new=c("tran", "gomf"))
#	tran2gobp  <- BioMart2df.fun(biomart="ensembl",  
#			dataset=dataset,
#			col.old=c("ensembl_transcript_id", "go_biological_pathway_id"),
#			col.new=c("tran", "gobp"))
#	tran2gocc  <- BioMart2df.fun(biomart="ensembl",  
#			dataset=dataset,
#			col.old=c("ensembl_transcript_id", "go_cellular_component_id"),
#			col.new=c("tran", "gocc"))
#}

# Get links between targets and GO terms from locally installed files:
if (summary.feature=="gene") {
	tran2gomf  <- read.table(file=file.path(path.mappings, "go_molecular_function_gene.txt"), sep="\t")
	tran2gobp  <- read.table(file=file.path(path.mappings, "go_biological_process_gene.txt"), sep="\t")
	tran2gocc  <- read.table(file=file.path(path.mappings, "go_cellular_component_gene.txt"), sep="\t")
}
if (summary.feature=="transcript") {
	tran2gomf  <- read.table(file=file.path(path.mappings, "go_molecular_function_transcript.txt"), sep="\t")
	tran2gobp  <- read.table(file=file.path(path.mappings, "go_biological_process_transcript.txt"), sep="\t")
	tran2gocc  <- read.table(file=file.path(path.mappings, "go_cellular_component_transcript.txt"), sep="\t")
}

# Fetch links between GO id and term from NCBI
#go2term <- GO2df.fun(url="ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz")

# Fetch links between GO id and term from locally installed files
go2term <- read.table(file="go_id_to_go_term.txt", sep="\t")

# Define the reference list to compare against, i.e. all transcripts in tran2gomf
if (summary.feature=="gene") {
	reference.list.mf <- as.vector(unique(tran2gomf[, "gene"]))
}
if (summary.feature=="transcript") {
	reference.list.mf <- as.vector(unique(tran2gomf[, "tran"]))
}
if (summary.feature=="gene") {
	reference.list.bp <- as.vector(unique(tran2gobp[, "gene"]))
}
if (summary.feature=="transcript") {
	reference.list.bp <- as.vector(unique(tran2gobp[, "tran"]))
}
if (summary.feature=="gene") {
	reference.list.cc <- as.vector(unique(tran2gocc[, "gene"]))
}
if (summary.feature=="transcript") {
	reference.list.cc <- as.vector(unique(tran2gocc[, "tran"]))
}

# Define the sample set, i.e. the set of targets in both targets.list and reference.list
sample.mf <- intersect(targets.list[,2], reference.list.mf)
sample.bp <- intersect(targets.list[,2], reference.list.bp)
sample.cc <- intersect(targets.list[,2], reference.list.cc)

# Run the hypergeometric test for over-representation
test.mf <- corna.test.fun(sample.mf,
		reference.list.mf,
		tran2gomf,
		fisher=T,
		sort="fisher",
		min.pop=5, 
		desc=go2term) 
test.bp <- corna.test.fun(sample.bp,
		reference.list.bp,
		tran2gobp,
		fisher=T,
		sort="fisher",
		min.pop=5, 
		desc=go2term) 
test.cc <- corna.test.fun(sample.cc,
		reference.list.cc,
		tran2gocc,
		fisher=T,
		sort="fisher",
		min.pop=5, 
		desc=go2term) 

# Fetch the list of significant terms
significant.mf <- test.mf[test.mf$fisher <= p.value.threshold, ]
significant.bp <- test.bp[test.bp$fisher <= p.value.threshold, ]
significant.cc <- test.cc[test.cc$fisher <= p.value.threshold, ]

# Write the results in three different tables
write.table(significant.mf, file="hyperg_go_mf.tsv", sep="\t")
write.table(significant.bp, file="hyperg_go_bp.tsv", sep="\t")
write.table(significant.cc, file="hyperg_go_cc.tsv", sep="\t")

