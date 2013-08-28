# TOOL ngs-correlation-analysis-mirna.R: "Correlate miRNA-seq and gene expression array data" (Performs a statistical test to detect miRNA targets whose expression is significantly positively or negatively correlated to the expression of the miRNA.)
# INPUT normalized_mirna.tsv: normalized_mirna.tsv TYPE GENE_EXPRS 
# INPUT normalized_gene.tsv: normalized_gene.tsv TYPE GENE_EXPRS 
# INPUT phenodata_mirna.tsv: phenodata_mirna.tsv TYPE GENERIC 
# INPUT phenodata_gene.tsv: phenodata_gene.tsv TYPE GENERIC 
# OUTPUT mirna-gene-positive-correlation.tsv: mirna-gene-positive-correlation.tsv 
# OUTPUT mirna-gene-negative-correlation.tsv: mirna-gene-negative-correlation.tsv 
# PARAMETER order.column.mirna: order.column.mirna TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing the order of the samples, so that the gene expression and miRNA expression arrays can be correctly matched in the analysis. For time course experiments the actual time can be used, for multiple-condition type of experiments it is adviced to encode the different conditions with a number, e.g. 1, 2, 3, 4 and 5 for an experiment where five different conditions have been assessed. NOTE: If a custom array was used for assessing the gene expression it is crucial that ENTREZ gene ID or HUGO gene symbols have been specified as identifier when importing the data into CHIPSTER.)
# PARAMETER order.column.gene: order.column.gene TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing the order of the samples, so that the gene expression and miRNA expression arrays can be correctly matched in the analysis. For time course experiments the actual time can be used, for multiple-condition type of experiments it is adviced to encode the different conditions with a number, e.g. 1, 2, 3, 4 and 5 for an experiment where five different conditions have been assessed. NOTE: If a custom array was used for assessing the gene expression it is crucial that ENTREZ gene ID have been specified as identifier when importing the data into CHIPSTER.)
# PARAMETER id.type: id.type TYPE [probe_id: probe_id, entrez_id: entrez_id] DEFAULT probe_id (Defines the type of gene identifier to use. For supported array types from Affymetrix, Agilent or Illumina probe_id should be used, whereas for custom arrays entrez_id should be used.)
# PARAMETER correlation.method: correlation.method TYPE [pearson: pearson, spearman: spearman, kendall: kendall] DEFAULT pearson (Method for calculating the correlation. Peasron's method is parametric, whereas Spearman's correlation is a non-parametric rank-based method that is less sensitive to outliers. Kendall's method is suitable in those cases one is interested in the sign of the changes in expression between adjacent data points, rather than the magnitude.)
# PARAMETER p.value.threshold: p.value.threshold TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER p.value.adjustment.method: p.value.adjustment.method TYPE [none: none, Bonferroni: Bonferroni, Holm: Holm, Hochberg: Hochberg, BH: BH, BY: BY] DEFAULT BH (Multiple testing correction method)

# Correlation analysis of miRNA targets
# MG, 11.2.2010
# IS, 8.6.2010, bugfix
# IS, 28.7.2010, now allows additional samples in the two data sets, but only those ones with a pair are used in the analysis
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

# Get sample order for matching the datasets
#mirna.order <- mirna.phenodata[,grep(order.column.mirna, colnames(mirna.phenodata))]
#gene.order <- gene.phenodata[,grep(order.column.gene, colnames(gene.phenodata))]

# check for unambiguity of sample identifiers
if (nrow(mirna.phenodata)!=length(unique(mirna.phenodata[,order.column.mirna])))
	stop('CHIPSTER-NOTE: Unambigous sample identifiers: ', paste(mirna.phenodata[,order.column.mirna], collapse=', ')) 
if (nrow(gene.phenodata)!=length(unique(gene.phenodata[,order.column.gene])))
	stop('CHIPSTER-NOTE: Unambigous sample identifiers: ', paste(gene.phenodata[,order.column.gene], collapse=', ')) 

# pick those samples that do have a matching pair
common.samples <- intersect(mirna.phenodata[,order.column.mirna], gene.phenodata[,order.column.gene])
rownames(mirna.phenodata) <- mirna.phenodata[,order.column.mirna]
rownames(gene.phenodata) <- gene.phenodata[,order.column.gene]
mirna.phenodata$n <- 1:nrow(mirna.phenodata)
gene.phenodata$n <- 1:nrow(gene.phenodata)
mirna.order <- mirna.phenodata[common.samples, 'n']
gene.order <- gene.phenodata[common.samples, 'n']

# Read the chiptype that was used for the gene expression data
if (id.type=="probe_id") {
	chip.type <- as.character(gene.phenodata[1,grep("chiptype", names(gene.phenodata))])
	if (length(grep(".db", chip.type)) == 0 & length(grep("pmcdf", chip.type)) == 0) {
		chip.type <- paste(chip.type, ".db", sep="")
	}
}
if (id.type=="entrez_id") {
	chip.type <- "org.Hs.eg.db"
}

# Sanity checks to make sure the experiment have enough conditions
if(length(unique(mirna.order))==1 | length(unique(gene.order))==1) {
	stop("You need to have at least 2 conditions or time points to run this analysis!")
}

# Sanity checks to make sure that the mirna and gene expression data sets
# have the same number of conditions
if(length(unique(mirna.order))!=length(unique(gene.order))) {
	stop("You need to have the same number of conditions, or time points, in the two data sets!")
}

# Define number of conditions
number.conditions <- length(mirna.order)

# Arrange the columns in the two datset so they match
mirna.data.3 <- mirna.data.2[,order(mirna.order)]
gene.data.3 <- gene.data.2[,order(gene.order)]

# Create data set appropriate for correlation testing
mirna.data.4 <- data.frame(mirna=rownames(mirna.data.3), exprs=mirna.data.3[,1])
gene.data.4 <- data.frame(gene=rownames(gene.data.3), exprs=gene.data.3[,1])

# check that the gene list actually contain at least one miRNA target
try(merged.table <- read.mir(gene=gene.data.4, mirna=mirna.data.4,
				annotation=chip.type), silent=TRUE)

if (match("merged.table",ls(),nomatch=0)==0) {
	stop("There were no targets found in either TarBase or PicTar databases for the supplied list of miRNA:s in the gene list selected. Try again by selecting a longer list of genes!")
}
merged.table <- read.mir(gene=gene.data.4, mirna=mirna.data.4,
		annotation=chip.type, verbose=TRUE)

for (count in 2:number.conditions) {
	mirna.data.4 <- data.frame(mirna=rownames(mirna.data.3), exprs=mirna.data.3[,count])
	gene.data.4 <- data.frame(gene=rownames(gene.data.3), exprs=gene.data.3[,count])
	temp.table <- read.mir(gene=gene.data.4, mirna=mirna.data.4,
			annotation=chip.type, verbose=TRUE)
	temp.table
	merged.table <- cbind (merged.table, temp.table)
}

# Change the symbols that come from the read.mir() function into
# the ones that come from the org.Hs.eg.db package
all_genes <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(all_genes)
mapped_genes  <- as.list(all_genes[mapped_genes])
symbols_list <- unlist(mapped_genes[as.character(merged.table$gene_id) ])
merged.table$symbol <- symbols_list

# Extract the matching mirna and gene expression values into two vectors
mirna.expression <- merged.table[, grep("mirExpr", names(merged.table))]
gene.expression <- merged.table[, grep("geneExpr", names(merged.table))]

# Calculate the pearson correlation value for each mirna-gene pair
results.table <- data.frame(merged.table$mature_miRNA, merged.table$gene_id, merged.table$symbol, correlation.coefficient=NA, correlation.p.value=NA)
names (results.table) <- c("miRNA", "entrez_id", "symbol", "correlation_coefficient", "p-value")
number.mirna <- dim(merged.table)[1]
for (mirna.count in 1:number.mirna) {
	correlation.coefficient <- cor (as.numeric(mirna.expression[mirna.count,]),as.numeric(gene.expression[mirna.count,]), method=correlation.method)
	correlation.p.value <- cor.test (as.numeric(mirna.expression[mirna.count,]),as.numeric(gene.expression[mirna.count,]), method=correlation.method)
	correlation.p.value <- correlation.p.value$p.value
	results.table[mirna.count,4] <- correlation.coefficient
	results.table[mirna.count,5] <- correlation.p.value
}
results.table[,5] <- p.adjust(results.table[,5], method=p.value.adjustment.method)

# Add rownames to allow use of Venn diagrams
# Construct the names from combining the miRNA nameand the target gene symbol
rownames (results.table) <- paste (results.table[,1],"_",results.table[,3],sep="")

# Find genes with statistically significant positive correlation
results.positive <- results.table[results.table[,4]>0,]
results.positive.significant <- results.positive[results.positive[,5]<=p.value.threshold,]

# Find genes with statistically significant negative correlation
results.negative <- results.table[results.table[,4]<=0,]
results.negative.significant <- results.negative[results.negative[,5]<=p.value.threshold,]

# Write the results to tables to be read into Chipster
write.table(results.positive.significant, file="mirna-gene-positive-correlation.tsv", sep="\t", quote=FALSE, row.names=TRUE)
write.table(results.negative.significant, file="mirna-gene-negative-correlation.tsv", sep="\t", quote=FALSE, row.names=TRUE)

# EOF
