# ANALYSIS Utilities/"Extract data for miRNA targets" (This tools provides a means of extracting the expression data
# for a set of target genes identified from an integrative analysis of miRNA and gene expression. Two files need to be selected,the
# result file from an integrative analysis and the file containing the gene expression data to be extracted.)
# INPUT GENE_EXPRS normalized_gene.tsv, GENERIC normalized_mirna.tsv, GENERIC phenodata_gene.tsv
# OUTPUT target-gene-expression.tsv, phenodata-target-gene.tsv
# PARAMETER common.column STRING DEFAULT empty (The name of the column that is common to the data tables.)

# Setting up test parameter settings
# common.column <- "symbol"

# Loads the normalized data and phenodata files
data_1 <- read.table(file="normalized_mirna.tsv", header=T, sep="\t", row.names=1)
data_2 <- read.table(file="normalized_gene.tsv", header=T, sep="\t", row.names=1)
phenodata <- read.table("phenodata_gene.tsv", header=T, sep="\t")

# Figure out which is the miRNA data
if (length (grep ("miRNA", names (data_1))) == 1) {
	mirna.data <- data_1
	gene.data <- data_2
}
if (length (grep ("miRNA", names (data_2))) == 1) {
	mirna.data <- data_2
	gene.data <- data_1
}

# Extract the common columns from the input data and condense
# to unique values to get rid of multiple entries
attach(mirna.data)
mirna.symbols <- get(common.column)
mirna.symbols <- unique(mirna.symbols)
detach(mirna.data)
attach(gene.data)
gene.symbols <- get(common.column)
detach(gene.data)
gene.data.2 <- gene.data [match(mirna.symbols, gene.symbols, nomatch=0),]

# write files
write.table(gene.data.2, file='target-gene-expression.tsv', quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
write.table(phenodata, file='phenodata-target-gene.tsv', quote=FALSE, sep='\t', na='', row.names=FALSE, col.names=TRUE)

# EOF

