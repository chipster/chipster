# ANALYSIS Pathways/"GO enrichment for miRNA targets" (Performs a statistical test for enrichments of GO terms in the predicted gene targets of a list of miRNA ID:s.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT hyperg_go.tsv
# PARAMETER ontology [all, biological_process, molecular_function, cellular_component] DEFAULT biological_process (the ontology to be analyzed)
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value threshold)
# PARAMETER p.adjust.method [none, BH, BY] DEFAULT BH (method for adjusting the p-value in order to account for multiple testing)
# PARAMETER over.or.under.representation [over, under] DEFAULT over (Should over or under-represented classes be seeked?)
# PARAMETER minimum.population INTEGER FROM 1 TO 1000000 DEFAULT 5 (minimum number of genes in in the reference list that map to a pathway)
# PARAMETER species [human, mouse, rat] DEFAULT human (the species for which the miRNA:s have been analyzed)

# POSSIBLE summary.feature [gene, transcript] DEFAULT gene (should the targets for the miRNA:s be transcripts or genes?)


# miRNA hypergeometric test for GO
# MG, 4.11.2009

# force "transcript" mode, since "gene" mode is not currently nicely handled
summary.feature <- "transcript"

# Reads the data
dat<-read.table("normalized.tsv", sep="\t", header=T)

# Extracts the identifiers
id<-as.character(rownames(dat))

# Translate parameter settings for biomaRt queries
if (species=="human") {
        dataset <- "hsapiens_gene_ensembl"
}
if (species=="mouse") {
        dataset <- "mmusculus_gene_ensembl"
}
if (species=="rat") {
        dataset <- "rnorvegicus_gene_ensembl"
}

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
        targets <- read.table(file.path(path.mappings, "mirna_mappings_rnorvegicus.txt"), sep="\t")
}

# Fetch the predicted targets for the list of miRNA:s
targets.list <- corna.map.fun(targets, id, "mir", "tran", all=TRUE)

# Get links between targets and GO terms from BioMart:
# Currently disabled to avoid connection problems to BiomaRT database
#if (summary.feature=="gene") {
#       tran2gomf  <- BioMart2df.fun(biomart="ensembl",  
#                       dataset=dataset,
#                       col.old=c("ensembl_gene_id", "go_molecular_function_id"),
#                       col.new=c("gene", "gomf"))
#       tran2gobp  <- BioMart2df.fun(biomart="ensembl",  
#                       dataset=dataset,
#                       col.old=c("ensembl_gene_id", "go_biological_process_id"),
#                       col.new=c("gene", "gobp"))
#       tran2gocc  <- BioMart2df.fun(biomart="ensembl",  
#                       dataset=dataset,
#                       col.old=c("ensembl_gene_id", "go_cellular_component_id"),
#                       col.new=c("gene", "gocc"))
#}
#if (summary.feature=="transcript") {
#       tran2gomf  <- BioMart2df.fun(biomart="ensembl",  
#                      dataset=dataset,
#                       col.old=c("ensembl_transcript_id", "go_molecular_function_id"),
#                       col.new=c("tran", "gomf"))
#       tran2gobp  <- BioMart2df.fun(biomart="ensembl",  
#                       dataset=dataset,
#                       col.old=c("ensembl_transcript_id", "go_biological_process_id"),
#                       col.new=c("tran", "gobp"))
#       tran2gocc  <- BioMart2df.fun(biomart="ensembl",  
#                       dataset=dataset,
#                       col.old=c("ensembl_transcript_id", "go_cellular_component_id"),
#                       col.new=c("tran", "gocc"))
#}

# Get links between targets and GO terms from locally installed files:
# Can be enabled if there are problems with BiomaRt interface
if (summary.feature=="gene" & species=="human") {
        tran2gomf  <- read.table(file=file.path(path.mappings, "go_molecular_function_gene_hsapiens.txt"), sep="\t")
        tran2gobp  <- read.table(file=file.path(path.mappings, "go_biological_process_gene_hsapiens.txt"), sep="\t")
        tran2gocc  <- read.table(file=file.path(path.mappings, "go_cellular_component_gene_hsapiens.txt"), sep="\t")
}
if (summary.feature=="transcript" & species=="human") {
        tran2gomf  <- read.table(file=file.path(path.mappings, "go_molecular_function_transcript_hsapiens.txt"), sep="\t")
        tran2gobp  <- read.table(file=file.path(path.mappings, "go_biological_process_transcript_hsapiens.txt"), sep="\t")
        tran2gocc  <- read.table(file=file.path(path.mappings, "go_cellular_component_transcript_hsapiens.txt"), sep="\t")
}
if (summary.feature=="gene" & species=="mouse") {
        tran2gomf  <- read.table(file=file.path(path.mappings, "go_molecular_function_gene_mmusculus.txt"), sep="\t")
        tran2gobp  <- read.table(file=file.path(path.mappings, "go_biological_process_gene_mmusculus.txt"), sep="\t")
        tran2gocc  <- read.table(file=file.path(path.mappings, "go_cellular_component_gene_mmusculus.txt"), sep="\t")
}
if (summary.feature=="transcript" & species=="mouse") {
        tran2gomf  <- read.table(file=file.path(path.mappings, "go_molecular_function_transcript_mmusculus.txt"), sep="\t")
        tran2gobp  <- read.table(file=file.path(path.mappings, "go_biological_process_transcript_mmusculus.txt"), sep="\t")
        tran2gocc  <- read.table(file=file.path(path.mappings, "go_cellular_component_transcript_mmusculus.txt"), sep="\t")
}
if (summary.feature=="gene" & species=="rat") {
        tran2gomf  <- read.table(file=file.path(path.mappings, "go_molecular_function_gene_rnorvegicus.txt"), sep="\t")
        tran2gobp  <- read.table(file=file.path(path.mappings, "go_biological_process_gene_rnorvegicus.txt"), sep="\t")
        tran2gocc  <- read.table(file=file.path(path.mappings, "go_cellular_component_gene_rnorvegicus.txt"), sep="\t")
}
if (summary.feature=="transcript" & species=="rat") {
        tran2gomf  <- read.table(file=file.path(path.mappings, "go_molecular_function_transcript_rnorvegicus.txt"), sep="\t")
        tran2gobp  <- read.table(file=file.path(path.mappings, "go_biological_process_transcript_rnorvegicus.txt"), sep="\t")
        tran2gocc  <- read.table(file=file.path(path.mappings, "go_cellular_component_transcript_rnorvegicus.txt"), sep="\t")
}

# Fetch links between GO id and term from NCBI, disabled for now due to unstable connections to NCBI
# Currently disabled to avoid connection problems to NCBI database
#go2term <- GO2df.fun(url="ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz")

# Fetch links between GO id and term from locally installed files
go2term <- read.table(file=file.path(path.mappings, "go_id_to_go_term.txt"), sep="\t")

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

# define the output variable
output <- data.frame(total=integer(0), expectation=numeric(0), observation=integer(0), p.value=numeric(0), description=character(0), ontology=character(0))

# test for over or under-represented categories
if (over.or.under.representation=='under') {
  hyper.lower.tail <- TRUE
} else {
  hyper.lower.tail <- FALSE
}

# Run the hypergeometric tests
if (ontology=="biological_process" || ontology=='all') {
  test.bp <- corna.test.fun(sample.bp,
                reference.list.bp,
                tran2gobp,
                hyper.lower.tail=hyper.lower.tail,
                fisher=F,
                sort="hypergeometric",
                min.pop=minimum.population,
                p.adjust.method=p.adjust.method,
                desc=go2term) 
  significant.bp <- test.bp[test.bp$hypergeometric <= p.value.threshold, ]
  if(nrow(significant.bp)>0) {
    significant.bp$ontology <- 'biological process'
    colnames(significant.bp) <- colnames(output)
    output <- rbind(output, significant.bp)
  }
}

if (ontology=="molecular_function" || ontology=='all') {
  test.mf <- corna.test.fun(sample.mf,
                reference.list.mf,
                tran2gomf,
                hyper.lower.tail=hyper.lower.tail,
                fisher=F,
                sort="hypergeometric",
                min.pop=minimum.population, 
                p.adjust.method=p.adjust.method,
                desc=go2term) 
  significant.mf <- test.mf[test.mf$hypergeometric <= p.value.threshold, ]
  if(nrow(significant.mf)>0) {
    significant.mf$ontology <- 'molecular function'
    colnames(significant.mf) <- colnames(output)
    output <- rbind(output, significant.mf)
  }
}

if (ontology=="cellular_component" || ontology=='all') {
  test.cc <- corna.test.fun(sample.cc,
                reference.list.cc,
                tran2gocc,
                hyper.lower.tail=hyper.lower.tail,
                fisher=F,
                sort="hypergeometric",
                min.pop=minimum.population,
                p.adjust.method=p.adjust.method,
                desc=go2term) 
  significant.cc <- test.cc[test.cc$hypergeometric <= p.value.threshold, ]
  if(nrow(significant.cc)>0) {
    significant.cc$ontology <- 'cellular component'
    colnames(significant.cc) <- colnames(output)
    output <- rbind(output, significant.cc)
  }
}

# write output
write.table(output, file="hyperg_go.tsv", sep="\t", quote=FALSE)

# EOF