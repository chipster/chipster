# TOOL ngs-up-down-analysis-mirna.R: "Up-down analysis of miRNA target genes with array data" (Given a dataset of miRNA expression and gene expression array data for a two-group comparison experiment, this tool identifies the genes whose expression is upregulated in response to an downregulated miRNA, or vice-versa. It is recommended that the two data sets have been filtered to exclude low quality and or invariable data and subjected to statistical testing for significant differences in expression between the two experiment groups. NOTE you need to assign a higher number to the samples belonging to the treatment grop in the column describing the experiment group for the calculations to be correct.)
# INPUT normalized_mirna.tsv: normalized_mirna.tsv TYPE GENE_EXPRS 
# INPUT normalized_gene.tsv: normalized_gene.tsv TYPE GENE_EXPRS 
# INPUT phenodata_mirna.tsv: phenodata_mirna.tsv TYPE GENERIC 
# INPUT phenodata_gene.tsv: phenodata_gene.tsv TYPE GENERIC 
# OUTPUT mirna-up-gene-down.tsv: mirna-up-gene-down.tsv 
# OUTPUT mirna-down-gene-up.tsv: mirna-down-gene-up.tsv 
# PARAMETER average.method: average.method TYPE [mean: mean, median: median] DEFAULT median (The method to calculate the average of samples in each experiment group.)
# PARAMETER groups.column.mirna: groups.column.mirna TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the experiment groups of the samples in the miRNA expression dataset.)
# PARAMETER groups.column.gene: groups.column.gene TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the experiment groups of the samples in the gene expression dataset.)
# PARAMETER id.type: id.type TYPE [probe_id: probe_id, entrez_id: entrez_id] DEFAULT probe_id (Defines the type of gene identifier to use. For supported array types from Affymetrix, Agilent or Illumina probe_id should be used, whereas for custom arrays entrez_id should be used.)

# Up-down analysis of miRNA targets
# MG, 25.2.2010
# MG, 10.2.2011, now gets the gene symbols from the org.HS.eg.db package and adds rownames the output tables to allow use of Venn diagram

# Loads the libraries
library(RmiR)
library(org.Hs.eg.db)

# Loads the normalized data and phenodata files
data_1 <- read.table(file="normalized_mirna.tsv", header=T, sep="\t", row.names=1)
data_2 <- read.table(file="normalized_gene.tsv", header=T, sep="\t", row.names=1)
phenodata_1 <- read.table("phenodata_mirna.tsv", header=T, sep="\t")
phenodata_2 <- read.table("phenodata_gene.tsv", header=T, sep="\t")

# Figure out which is the miRNA data
if (phenodata_1$chiptype[1] == "miRNA") {
	mirna.phenodata <- phenodata_1
	mirna.data <- data_1
	gene.phenodata <- phenodata_2
	gene.data <- data_2
}
if (phenodata_2$chiptype[1] == "miRNA") {
	mirna.phenodata <- phenodata_2
	mirna.data <- data_2
	gene.phenodata <- phenodata_1
	gene.data <- data_1
}

# Separates expression values and flags
mirna.data.2 <- mirna.data[,grep("chip", names(mirna.data))]
gene.data.2 <- gene.data[,grep("chip", names(gene.data))]

# Get experiment groups from for the two datasets
mirna.groups <- mirna.phenodata[,grep(groups.column.mirna, colnames(mirna.phenodata))]
gene.groups <- gene.phenodata[,grep(groups.column.gene, colnames(gene.phenodata))]

# Sanity checks to make sure only two experiment groups exist in the two datasets
if(length(unique(mirna.groups))==1 | length(unique(mirna.groups))>=3) {
	stop("You need to have exactly two groups in the miRNA expression dataset to run this analysis")
}
if(length(unique(gene.groups))==1 | length(unique(gene.groups))>=3) {
	stop("You need to have exactly two groups in the gene expression dataset to run this analysis")
}

# Calculate the average expression for the different experiment groups
mirna.average.1 <- mirna.data.2[,mirna.groups==sort(unique(mirna.groups), decreasing=TRUE)[1]]
mirna.average.1 <- apply(mirna.average.1, FUN=average.method, MARGIN=1)
mirna.average.2 <- mirna.data.2[,mirna.groups==sort(unique(mirna.groups), decreasing=TRUE)[2]]
mirna.average.2 <- apply(mirna.average.2, FUN=average.method, MARGIN=1)
mirna.ratio <- mirna.average.1-mirna.average.2
gene.average.1 <- gene.data.2[,gene.groups==sort(unique(gene.groups), decreasing=TRUE)[1]]
gene.average.1 <- apply(gene.average.1, FUN=average.method, MARGIN=1)
gene.average.2 <- gene.data.2[,gene.groups==sort(unique(gene.groups),decreasing=TRUE)[2]]
gene.average.2 <- apply(gene.average.2, FUN=average.method, MARGIN=1)
gene.ratio <- gene.average.1-gene.average.2

# Read the chiptype that was used for the gene expression data
if (id.type=="probe_id") {
	chip.type <- as.character(gene.phenodata[1,grep("chiptype", names(gene.phenodata))])
}
if (id.type=="entrez_id") {
	chip.type <- "org.Hs.eg.db"
}

# Construct datasets suitable for read.mir() function
mirna.data.3 <- cbind(names(mirna.ratio),as.numeric(mirna.ratio))
gene.data.3 <- cbind(names(gene.ratio),as.numeric(gene.ratio))
mirna.data.4 <- as.data.frame(mirna.data.3)
gene.data.4 <- as.data.frame(gene.data.3)
mirna.data.4[,2] <- as.numeric(mirna.data.3[,2])
gene.data.4[,2] <- as.numeric(gene.data.3[,2])

# Check that the gene list actually contain at least one miRNA target
try(merged.table <- read.mir(gene=gene.data.4, mirna=mirna.data.4,
				annotation=chip.type), silent=FALSE)
if (match("merged.table",ls(),nomatch=0)==0) {
	stop("There were no targets found in either TarBase or PicTar databases
					for the supplied list of miRNA:s in the gene list selected.
					Try again by selecting a longer list of genes!")
}
merged.table <- read.mir(gene=gene.data.4, mirna=mirna.data.4,
		annotation=chip.type, verbose=TRUE)

# Change the symbols that come from the read.mir() function into
# the ones that come from the org.Hs.eg.db package
all_genes <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(all_genes)
mapped_genes  <- as.list(all_genes[mapped_genes])
symbols_list <- unlist(mapped_genes[as.character(merged.table$gene_id) ])
merged.table$symbol <- symbols_list

# Change the column names for sake of consistency with other tools
colnames (merged.table) [1] <- "entrez_id"
colnames (merged.table) [2] <- "miRNA"

# Add rownames to allow use of Venn diagrams
# Construct the names from combining the miRNA nameand the target gene symbol
rownames (merged.table) <- paste (merged.table[,2],"_",merged.table[,4],sep="")

# Identify the miRNA-gene pairs that are behaving oppositely
mirna.ratio <- merged.table$mirExpr
gene.ratio <- merged.table$geneExpr
up.mirna.down.gene <- merged.table[(mirna.ratio>=0 & gene.ratio<0),1:5]
down.mirna.up.gene <- merged.table[(mirna.ratio<0 & gene.ratio>=0),1:5]

# Write the results to tables to be read into Chipster
write.table(up.mirna.down.gene, file="mirna-up-gene-down.tsv", sep="\t", quote=FALSE, row.names=TRUE)
write.table(down.mirna.up.gene, file="mirna-down-gene-up.tsv", sep="\t", quote=FALSE, row.names=TRUE)

# EOF
