# TOOL dea-deseq-RNA.R: "Differential expression analysis using DESeq" (This tool will perform an analysis for differentially expressed sequences using the R implementation of the DESeq algorithm.)
# INPUT data.tsv TYPE GENERIC
# INPUT phenodata.tsv TYPE GENERIC
# OUTPUT OPTIONAL de-list.tsv
# OUTPUT OPTIONAL de-list.bed
# OUTPUT OPTIONAL ma-plot.pdf
# OUTPUT OPTIONAL dispersion-plot.pdf
# OUTPUT OPTIONAL p-value-plot.pdf
# PARAMETER column: "Column describing groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER normalization: "Apply normalization" TYPE [yes, no] DEFAULT yes (If enabled, a normalization factor based on estimated library size is calculated.)
# PARAMETER dispersion_method: "Dispersion method" TYPE [common, tagwise] DEFAULT tagwise (The dispersion of counts for any given sequence can either be estimated based on the actual counts in the sample data set or be moderated across a selection of sequences with similar count numbers. The latter option, which is set by default, typically yields higher sensitivity and specificity. Note that when no biological replicates are available common dispersion is used regardless of the setting.)
# PARAMETER dispersion_estimate:"Dispersion estimate" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.1 (The value to use for estimating the common dispersion when no replicates are available.) 
# PARAMETER p.value.adjustment.method: "Multiple testing correction" TYPE [none, Bonferroni, Holm, Hochberg, BH, BY] DEFAULT BH (Multiple testing correction method.)
# PARAMETER p.value.cutoff: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (The cutoff for statistical significance.)
# PARAMETER image_width: "Plot width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image_height: "Plot height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


############################################################
#                                                          #
# Analaysis workflow using DESeq for normalization and     #
# statistical testing for finding differentially expressed #
# sequence tags                                            #
#                                                          #
# MG, 7.2.2012                                             #
#                                                          #
############################################################

# Loads the libraries
library(DESeq)

# Simplify variable names
w <- image_width
h <- image_height

# Set parameters for testing
# p.value.cutoff <- 0.1

# Loads the normalized data
file <- c("data.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
annotations <- dat[,-grep("chip", names(dat))]
dat2 <- dat[,grep("chip", names(dat))]

# Test needs a parameter "groups" that specifies the grouping of the samples
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
groups <- as.character (phenodata[,pmatch(column,colnames(phenodata))])
group_levels <- levels(as.factor(groups))
number_samples <- length(groups)

# If the library_size column contains data then use that as estimate
lib_size <- as.numeric(phenodata$library_size)
if (is.na(lib_size[1])) estimate_lib_size <- "TRUE" else estimate_lib_size <- "FALSE"

# Sanity checks
# only 2 group comparison is supported
if (length(unique(groups))==1 | length(unique(groups))>=3) {
	stop("CHIPSTER-NOTE: You need to have exactly two groups to run this analysis")
}
# if no biological replicates, force common dispersion
#if (number_samples == 2) dispersion_method <- "common" 

# Create a counts data object
counts_data <- newCountDataSet( dat2, groups )

# Calculate scaling factors based on estimated library size
counts_data <- estimateSizeFactors( counts_data )

# Estimate dispersion values for each gene and replaced with fitted values
# use sharingMode parameter to control how conservative
# use fitType to control for parametric or local fit
counts_data <- estimateDispersions( counts_data )

# Function that produces a qc plot to check dispersion estimates
plotDispEsts <- function( cds ) {
	plot(rowMeans( counts( cds, normalized=TRUE ) ), fitInfo(cds)$perGeneDispEsts,pch = '.', 
			log="xy", main="Dispersion plot", xlab="normalized counts", ylab="dispersion" )
	xg <- 10^seq( -.5, 5, length.out=300 )
	lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
	legend (x="topright", legend="fitted dipersion", col="red", cex=1, pch="-")
}

# Make plot
pdf(file="dispersion-plot.pdf")
plotDispEsts( counts_data )
dev.off()

# Calculate statistic for differential expression
results_table <- nbinomTest( counts_data, group_levels[2], group_levels[1] )

# Filter out the significant ones
significant_table <- results_table[ (results_table$padj <  p.value.cutoff ),]

# Order results based on raw p-values
significant_table <- significant_table[ order(significant_table$pval), ] 

# Output the table
if (dim(significant_table)[1] > 0) {
	write.table(significant_table, file="de-list.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# Make histogram of p-values
pdf (file="p-value-plot.pdf")
hist(results_table$pval, breaks=100, col="blue", border="slateblue", freq=FALSE,
		main="P-value distribution", xlab="p-value", ylab="proportion (%)")
hist(results_table$padj, breaks=100, col="red", border="slateblue", add=TRUE, freq=FALSE)
abline(h=1, lwd=2, lty=2, col="black")
abline(v=p.value.cutoff, lwd=2, lty=2, col="green")
legend (x="topright", legend=c("p-values","adjusted p-values", "uniform distribution", "significance cutoff"), col=c("blue","red","black","green"),
		cex=1, pch=15)
dev.off()



# Define function for making MA-plot of significant findings
plotDE <- function( res )
	plot(res$baseMean, res$log2FoldChange,
			log="x", pch=20, cex=.3,  col = ifelse( res$padj < .1, "red", "black"),
			main="MA plot", xlab="mean counts", ylab="log2(fold change)" ) 
# Make MA-plot
pdf(file="ma-plot.pdf")
plotDE( results_table )
legend (x="topleft", legend=c("significant features","not significant"), col=c("red","black"),
		cex=1, pch=19)
dev.off()
# EOF







# Create a tbale with the original counts per sample together with the statistical tests results
# ready for output in Chipster
# If there are no significant results return a message
if (dim(significant_results)[1] > 0) {
	output_table <- data.frame (dat[significant_indices,], significant_results)
}

# Output the table
if (dim(significant_results)[1] > 0) {
	write.table(output_table, file="de-list.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# Also output a bed graph file for visualization and region matching tools
if (dim(significant_results)[1] > 0) {
	empty_column <- character(length(significant_indices))
	bed_output <- output_table [,c("chr","start","end")]
	bed_output <- cbind(bed_output,empty_column)
	bed_output <- cbind(bed_output, output_table[,"logFC"])
	write.table(bed_output, file="de-list.bed", sep="\t", row.names=F, col.names=F, quote=F)
}

# Output a message if no significant genes are found
if (dim(significant_results)[1] == 0) {
	cat("No statistically significantly expressed sequences were found. Try again with a less stringent p-value cut-off or multiple testing correction method.", file="edgeR-log.txt")
}

# EOF