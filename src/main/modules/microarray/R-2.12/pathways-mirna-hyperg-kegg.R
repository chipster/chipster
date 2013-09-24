# TOOL pathways-mirna-hyperg-kegg.R: "KEGG enrichment for miRNA targets" (Given a list of miRNA identifiers, tests for enrichment of KEGG pathways in their predicted gene targets.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT hyperg_kegg.tsv: hyperg_kegg.tsv 
# PARAMETER p.value.threshold: p.value.threshold TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value threshold)
# PARAMETER p.adjust.method: p.adjust.method TYPE [none: none, BH: BH, BY: BY] DEFAULT BH (Method for adjusting the p-value in order to account for multiple testing)
# PARAMETER minimum.population: minimum.population TYPE INTEGER FROM 1 TO 1000000 DEFAULT 10 (Minimum number of genes in in the reference list that map to a pathway)
# PARAMETER species: species TYPE [human: human, mouse: mouse, rat: rat] DEFAULT human (The species for which the miRNA:s have been analyzed)

# POSSIBLE summary.feature [gene, transcript] DEFAULT gene (should the targets for the miRNA:s be transcripts or genes?)

# miRNA hypergeometric test for KEGG
# MG, 4.11.2009
# modifed 16.12.2009 by MG

# force "transcript" mode
summary.feature <- "transcript"

# Reads the data
dat<-read.table("normalized.tsv", sep="\t", header=T)

# Extracts the identifiers
id<-as.character(rownames(dat))

# Translate parameter settings for biomaRt queries
if (species=="human") {
	dataset <- "hsapiens_gene_ensembl"
}
#if (species=="human" & summary.feature=="transcript") {
#        dataset <- "hsapiens_transcript_ensembl"
#}
if (species=="mouse") {
	dataset <- "mmusculus_gene_ensembl"
}
#if (species=="mouse" & summary.feature=="transcript") {
#        dataset <- "mmusculus_t_ensembl"
#}
if (species=="rat") {
	dataset <- "rnorvegicus_gene_ensembl"
}
#if (species=="rat" & summary.feature=="transcript") {
#        dataset <- "rnorvegicus_transcript_ensembl"
#}


# Read in the CORNA library, which contains the functions to map miRNA:s to targets
# and performs hypergeometric test for enrichment of GO terms
library(CORNA)

# Download the mapping of miRNA to its targets from Sanger institute
# Currently disabled becuse of unstable web-services
#if (species=="human") {
#	targets <- miRBase2df.fun(url="ftp://ftp.sanger.ac.uk/pub/mirbase/targets/v5/arch.v5.txt.homo_sapiens.zip")
#}
#if (species=="mouse") {
#	targets <- miRBase2df.fun(url="ftp://ftp.sanger.ac.uk/pub/mirbase/targets/v5/arch.v5.txt.mus_musculus.zip")
#}
#if (species=="rat") {
#	targets <- miRBase2df.fun(url="ftp://ftp.sanger.ac.uk/pub/mirbase/targets/v5/arch.v5.txt.rattus_norvegicus.zip")
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


# obtain a link from transcript to gene from BiomaRt
# disabled for now to avoid connection problems to BiomaRt
#tran2gene <- BioMart2df.fun(
#		biomart="ensembl",
#		dataset=dataset,
#		col.old=c("ensembl_transcript_id",
#				"ensembl_gene_id"),
#		col.new=c("tran", "gene"))

# obtain a link from transcript to gene from locally installed files
if (species=="human") {
	tran2gene <- read.table(file=file.path(path.mappings,"transcripts_to_genes_hsapiens.txt"), sep="\t")
}
if (species=="mouse") {
	tran2gene <- read.table(file=file.path(path.mappings,"transcripts_to_genes_mmusculus.txt"), sep="\t")
}
if (species=="rat") {
	tran2gene <- read.table(file=file.path(path.mappings,"transcripts_to_genes_rnorvegicus.txt"), sep="\t")
}

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

# read pathway information from KEGG
# disabled for now to avoid connection problems to KEGG
#if (species="human") {
#	gene2path <- KEGG2df.fun(org="hsa")
#}
#if (species="mouse") {
#	gene2path <- KEGG2df.fun(org="mmu")
#}
#if (species="rat") {
#	gene2path <- KEGG2df.fun(org="rno")
#}
	
# read pathway information from local files
if (species=="human") {
	gene2path <- read.table(file=file.path(path.mappings,"genes_to_kegg_hsapiens.txt"), sep="\t")
}
if (species=="mouse") {
	gene2path <- read.table(file=file.path(path.mappings,"genes_to_kegg_mmusculus.txt"), sep="\t")
}
if (species=="rat") {
	gene2path <- read.table(file=file.path(path.mappings,"genes_to_kegg_rnorvegicus.txt"), sep="\t")
}
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
		fisher=F,
		chi.square=F,
		#fisher.alternative="greater",
		min.pop=minimum.population,
		sort="hypergeometric",
		p.adjust.method=p.adjust.method,
		desc=path2name)

# fetch significant pathways
significant.kegg <- test[test$hypergeometric<=p.value.threshold,]

# write results table
write.table(significant.kegg, file="hyperg_kegg.tsv", sep="\t", quote=F)


