# TOOL dea-deseq.R: "Differential expression using DESeq" (Differential expression analysis using the DESeq Bioconductor package. You can create the input count table and phenodata file using the tool \"Utilities - Define NGS experiment\".)
# INPUT data.tsv TYPE GENERIC
# INPUT META phenodata.tsv TYPE GENERIC
# OUTPUT OPTIONAL de-list-deseq.tsv
# OUTPUT OPTIONAL de-list-deseq.bed
# OUTPUT OPTIONAL ma-plot-deseq.pdf
# OUTPUT OPTIONAL dispersion-plot-deseq.pdf
# OUTPUT OPTIONAL p-value-plot-deseq.pdf
# PARAMETER column: "Column describing groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test.)
# PARAMETER OPTIONAL normalization: "Apply normalization" TYPE [yes, no] DEFAULT yes (Should effective library size be estimated. This corrects for RNA composition bias. Note that if you have supplied library size in phenodata, size factors are calculated based on the library size total, and composition bias is not corrected.)
# PARAMETER OPTIONAL dispersion_estimate:"Dispersion estimation method" TYPE [parametric: "parametric", local: "local"] DEFAULT parametric (Dispersion can be estimated using a local fit or a two-coefficient parametric model. You should use local fit if there are no biological replicates.)
# PARAMETER OPTIONAL fitting_method: "Use fitted dispersion values" TYPE [maximum: "when higher than original values", fit-only: "always"] DEFAULT maximum (Should the dispersion of counts for a gene be replaced with the fitted value always, or only when the fitted value is larger? Replacing always optimises the balance between false positives and false negatives. Replacing only when the fitted value is higher is more conservative and minimizes false positives.)
# PARAMETER OPTIONAL p.value.adjustment.method: "Multiple testing correction" TYPE [none, bonferroni: "Bonferroni", holm: "Holm", hochberg: "Hochberg", BH: "BH", BY: "BY"] DEFAULT BH (Multiple testing correction method.)
# PARAMETER OPTIONAL p.value.cutoff: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (The cutoff for adjusted p-value.)
# PARAMETER OPTIONAL image_width: "Plot width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image.)
# PARAMETER OPTIONAL image_height: "Plot height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image.)

                                                          
# MG 7.2.2012                                             
# EK 6.5.2012, clarified texts
# EK 12.5.2012, fixed the fitting method parameter
# EK 30.4.2013, added BED sorting, made genomic location info optional so that external count tables can be used
# EK 6.5.2013, removed replicates parameter
# MK 29.01.2013, fixed bug why FoldChange column was duplicated in results

# Loads the libraries
library(DESeq)

# Loads the counts data
file <- c("data.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)

# Extracts expression values
annotations <- dat[,-grep("chip", names(dat))]
dat2 <- dat[,grep("chip", names(dat))]

# Test needs a parameter "groups" that specifies the grouping of the samples
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
groups <- as.character (phenodata[,pmatch(column,colnames(phenodata))])
group_levels <- levels(as.factor(groups))
number_samples <- length(groups)

# If the library_size column of phenodata contains data, then use that to estimate size factors
lib_size <- as.numeric(phenodata$library_size)
if (is.na(lib_size[1])) estimate_lib_size <- "TRUE" else estimate_lib_size <- "FALSE"
lib_size <- lib_size/mean(lib_size)

# Sanity check: only 2 group comparison is supported
if (length(unique(groups))==1 | length(unique(groups))>=3) {
	stop("CHIPSTER-NOTE: You need to have exactly two groups to run this analysis")
}

# Create a counts data object
counts_data <- newCountDataSet(dat2, groups)

# Calculate size factors based on estimated library size, unless it is given in phenodata
# Set size factors to 1 if normalization is turned off
if (normalization == "yes") {
	if (estimate_lib_size) {
		counts_data <- estimateSizeFactors(counts_data)
	} else {
		counts_data <- estimateSizeFactors(counts_data)
		sizeFactors(counts_data) <- lib_size
	}
} else {
	sizeFactors(counts_data) <- 1
}

# For later use: filter out genes which have less than 5 counts in user-defined number of samples 
#if (filter > 0) {
#	keep <- rowSums(counts(counts_data) >5) >= filter
#	counts_data <- counts_data[keep,]
#}

# Estimate dispersion values for each gene and replace with fitted values always or only when the fitted value is higher
# Use fitType to control for parametric or local fit
if (number_samples == 2 ) {
	counts_data <- estimateDispersions(counts_data, method="blind", sharingMode="fit-only", fitType=dispersion_estimate)
} else {
	counts_data <- estimateDispersions(counts_data, method="pooled", sharingMode=fitting_method, fitType=dispersion_estimate)
}

# Function that produces a dispersion plot
plotDispEsts <- function(cds) {
	plot(rowMeans( counts(cds, normalized=TRUE)), fitInfo(cds)$perGeneDispEsts,pch = '.', 
			log="xy", main="Dispersion plot", xlab="normalized counts", ylab="dispersion")
	xg <- 10^seq( -.5, 5, length.out=300)
	lines(xg, fitInfo(cds)$dispFun(xg), col="red")
	legend(x="topright", legend="fitted dispersion", col="red", cex=1, pch="-")
}

# Make dispersion plot
pdf(file="dispersion-plot-deseq.pdf")
plotDispEsts(counts_data)
dev.off()

# Calculate statistic for differential expression
results_table <- nbinomTest(counts_data, group_levels[1], group_levels[2] )

# Merge with original data table

output_table <- cbind (dat, results_table[,-1])

# Adjust p-values
output_table$padj <- p.adjust(output_table$pval, method=p.value.adjustment.method)

# Keep significant DEGs
significant_table <- output_table[ (output_table$padj <  p.value.cutoff),]

# Remove rows with NA adjusted p-values
significant_table <- significant_table[! (is.na(significant_table$padj)),]

# Order results based on raw p-values
significant_table <- significant_table[ order(significant_table$pval), ] 


#############################


#significant_indices <- rownames (significant_table)

# If significant results are found, create an output table with the original counts per sample together with the statistical tests results
# If genomic coordinates are present, output a sorted BED file for genome browser visualization and region matching tools
#if (dim(significant_table)[1] > 0) {
#	output_table <- data.frame (dat[significant_indices,], significant_table)
#	write.table(output_table, file="de-list-deseq.tsv", sep="\t", row.names=T, col.names=T, quote=F)
#	these.colnames <- colnames(dat)
#	if("chr" %in% these.colnames) {
#		empty_column <- character(length(significant_indices))
#		bed_output <- output_table [,c("chr","start","end")]
#		bed_output <- cbind(bed_output,empty_column)
#		bed_output <- cbind(bed_output, output_table[,"log2FoldChange"])
#		source(file.path(chipster.common.path, "bed-utils.R"))
#		bed_output <- sort.bed(bed_output)
#		write.table(bed_output, file="de-list-deseq.bed", sep="\t", row.names=F, col.names=F, quote=F)
#	}	
#}

#############################


# Output the table
if (dim(significant_table)[1] > 0) {
	ndat <- ncol(dat)
	nmax <- ncol(significant_table)
	write.table(cbind(significant_table[,1:ndat], round(significant_table[, (ndat+1):(nmax-2)], digits=2), format(significant_table[, (nmax-1):nmax], digits=4, scientific=T)), file="de-list-deseq.tsv", sep="\t", row.names=T, col.names=T, quote=F)

	#write.table(significant_table, file="de-list-deseq.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# If genomic coordinates are present, output a sorted BED file for genome browser visualization and region matching tools
source(file.path(chipster.common.path, "bed-utils.R"))
these.colnames <- colnames(dat)
if("chr" %in% these.colnames) {
	if (dim(significant_table)[1] > 0) {
		empty_column <- character(length(significant_table[1]))
		bed_output <- significant_table [,c("chr","start","end")]
		bed_output <- cbind(bed_output,empty_column)
		bed_output <- cbind(bed_output, significant_table[,"log2FoldChange"])
		bed_output <- sort.bed(bed_output)
		write.table(bed_output, file="de-list-deseq.bed", sep="\t", row.names=F, col.names=F, quote=F)
	}
}

 
# Make histogram of p-values with overlaid significance cutoff and uniform distribution
pdf (file="p-value-plot-deseq.pdf")
hist(output_table$pval, breaks=100, col="blue",
		border="slateblue", freq=FALSE,
		main="P-value distribution", xlab="p-value", ylab="proportion (%)")
hist(output_table$padj, breaks=100, col="red",
		border="slateblue", add=TRUE, freq=FALSE)
abline(h=1, lwd=2, lty=2, col="black")
abline(v=p.value.cutoff, lwd=2, lty=2, col="green")
legend (x="topright", legend=c("p-values","adjusted p-values", "uniform distribution", "significance cutoff"), col=c("blue","red","black","green"),
		cex=1, pch=15)
dev.off()

# Define function for making MA-plot of significant findings
plotDE <- function(res)
	plot(res$baseMean, res$log2FoldChange,
			log="x", pch=20, cex=.25, col = ifelse( res$padj < p.value.cutoff, "red", "black"),
			main="MA plot", xlab="mean counts", ylab="log2(fold change)") 

# Make MA-plot
pdf(file="ma-plot-deseq.pdf")
plotDE(output_table)
legend (x="topleft", legend=c("significant","not significant"), col=c("red","black"),
		cex=1, pch=19)
abline(h = c(-1, 0, 1), col = c("dodgerblue", "darkgreen", "dodgerblue"), lwd = 2)
dev.off()

# EOF
