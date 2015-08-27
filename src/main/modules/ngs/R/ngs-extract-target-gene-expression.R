# TOOL ngs-extract-target-gene-expression.R: "Extract gene expression array data for miRNA targets" (Extracts expression array data for genes which correlate with miRNA expression, as judged by an integrative analysis of miRNA and gene expression. Two files need to be selected: the result file from the correlation analysis and the file containing the gene expression array data.)
# INPUT normalized_gene.tsv: normalized_gene.tsv TYPE GENE_EXPRS 
# INPUT normalized_mirna.tsv: normalized_mirna.tsv TYPE GENERIC 
# INPUT META phenodata_gene.tsv: phenodata_gene.tsv TYPE GENERIC 
# OUTPUT target-gene-expression.tsv: target-gene-expression.tsv 
# OUTPUT META phenodata-target-gene.tsv: phenodata-target-gene.tsv 
# PARAMETER common.column: "The gene name column in both files" TYPE STRING DEFAULT empty (The gene name column that is common to both data tables.)

# MK 15.05.2014, modified to R3

# Loads the normalized data and phenodata files
data_1 <- read.table(file="normalized_mirna.tsv", header=T, sep="\t", row.names=1)
data_2 <- read.table(file="normalized_gene.tsv", header=T, sep="\t", row.names=1)
phenodata <- read.table("phenodata_gene.tsv", header=T, sep="\t")

# Figure out which is the miRNA data
if (length (grep ("mirna", names (data_1), ignore.case = TRUE)) == 1) {
	mirna.data <- data_1
	gene.data <- data_2
}
if (length (grep ("mirna", names (data_2), ignore.case = TRUE)) == 1) {
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

