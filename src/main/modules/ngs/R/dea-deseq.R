# TOOL dea-deseq.R: "Differential expression using DESeq" (Differential expression analysis using the DESeq Bioconductor package. You can create the input count table and phenodata file using the tool \"Utilities - Define NGS experiment\".)
# INPUT data.tsv TYPE GENERIC
# INPUT META phenodata.tsv TYPE GENERIC
# OUTPUT OPTIONAL de-list-deseq.tsv
# OUTPUT OPTIONAL de-list-deseq.bed
# OUTPUT OPTIONAL deseq_report.pdf
# PARAMETER column: "Column describing groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test.)
# PARAMETER OPTIONAL ad_factor: "Column describing additional experimental factor" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing an additional experimental factor. If given, p-values in the output table are from a likelihood ratio test of a model including the experimental groups and experimental factor vs a model which only includes the experimental factor.)
# PARAMETER OPTIONAL normalization: "Apply normalization" TYPE [yes, no] DEFAULT yes (Should effective library size be estimated. This corrects for RNA composition bias. Note that if you have supplied library size in phenodata, size factors are calculated based on the library size total, and composition bias is not corrected.)
# PARAMETER OPTIONAL dispersion_estimate:"Dispersion estimation method" TYPE [parametric: "parametric", local: "local"] DEFAULT parametric (Dispersion can be estimated using a local fit or a two-coefficient parametric model. You should use local fit if there are no biological replicates.)
# PARAMETER OPTIONAL fitting_method: "Use fitted dispersion values" TYPE [maximum: "when higher than original values", fit-only: "always", gene-est-only: "no fitting"] DEFAULT maximum (Should the dispersion of counts for a gene be replaced with the fitted value always, or only when the fitted value is larger? Replacing always optimises the balance between false positives and false negatives. Replacing only when the fitted value is higher is more conservative and minimizes false positives. No fitting is preferable when the number of replicates is large.)
# PARAMETER OPTIONAL p.value.adjustment.method: "Multiple testing correction" TYPE [none, bonferroni: "Bonferroni", holm: "Holm", hochberg: "Hochberg", BH: "BH", BY: "BY"] DEFAULT BH (Multiple testing correction method.)
# PARAMETER OPTIONAL p.value.cutoff: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (The cutoff for adjusted p-value.)

                                                          
# MG 7.2.2012                                             
# EK 6.5.2012, clarified texts
# EK 12.5.2012, fixed the fitting method parameter
# EK 30.4.2013, added BED sorting, made genomic location info optional so that external count tables can be used
# EK 6.5.2013, removed replicates parameter
# MK 29.01.2013, fixed bug why FoldChange column was duplicated in results
# MK 14.04.2014, added possibility to compare GLM models
# MK 15.04.2014, added possibility to use DESeq2
# AMS 17.06.2014, split DESeq2 to a separate tool
# EK 30.6.2014, added gene names to BED output, clarified the script
# AMS 07.04.2015, joined PDF outputs

# Loads the libraries
source(file.path(chipster.common.path, "bed-utils.R"))
library(DESeq)

# Loads the counts data and extracts expression values
file <- c("data.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)
dat2 <- dat[,grep("chip", names(dat))]

# Test needs a parameter "groups" that specifies experimental conditions
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
groups <- as.character (phenodata[,pmatch(column,colnames(phenodata))])
group_levels <- levels(as.factor(groups))

# Sanity check: only 2 group comparison is supported
if (length(unique(groups)) == 1 | (length(unique(groups)) >= 3)) {
	stop("CHIPSTER-NOTE: You need to have exactly two groups to run DESeq. Use DESeq2 if you have multifactor, multigroup data")
}

# Read the additional factor for GLM test. Construct design matrix from this 
exp_factor <- NULL
if(ad_factor != "EMPTY") {
	exp_factor <- as.character (phenodata[,pmatch(ad_factor,colnames(phenodata))])
	design <- data.frame(condition=as.factor(groups), exp_factor=exp_factor)
	rownames(design) <- colnames(dat2)
}

# Create a counts data object
if(ad_factor == "EMPTY") {
	counts_data <- newCountDataSet(dat2, groups)
} else if(ad_factor != "EMPTY") {
	counts_data <- newCountDataSet(dat2, design)
}

# Calculate size factors based on estimated library size, unless library sizes are given in phenodata
# Set size factors to 1 if normalization is turned off
lib_size <- as.numeric(phenodata$library_size)
if (normalization == "yes" && is.na(lib_size[1])) {
	counts_data <- estimateSizeFactors(counts_data)
} else if(normalization == "yes" && !(is.na(lib_size[1]))) {
	lib_size <- lib_size/mean(lib_size)
	counts_data <- estimateSizeFactors(counts_data)
	sizeFactors(counts_data) <- lib_size
} else {
	sizeFactors(counts_data) <- 1
}

# Estimate dispersion values for each gene and replace with fitted values always or only when the fitted value is higher
# Use fitType to control for parametric or local fit. Method blind can be used when data has no biological replicates. 
if(ad_factor == "EMPTY" && length(groups) == 2) {
	counts_data <- estimateDispersions(counts_data, method="blind", sharingMode="fit-only", fitType=dispersion_estimate)
} else if(ad_factor == "EMPTY" && length(groups) > 2) {
	counts_data <- estimateDispersions(counts_data, method="pooled", sharingMode=fitting_method, fitType=dispersion_estimate)
} else if(ad_factor != "EMPTY") {
	counts_data <- estimateDispersions(counts_data, method="pooled-CR", sharingMode=fitting_method, fitType=dispersion_estimate, modelFormula = count ~ exp_factor + condition)	
}

# Calculate statistic for differential expression
if(ad_factor == "EMPTY") {
	results_table <- nbinomTest(counts_data, group_levels[1], group_levels[2] )[,-1]
} else if (ad_factor != "EMPTY") {
	results_table <- nbinomTest(counts_data, group_levels[1], group_levels[2] )[,-1]

	fit0 = fitNbinomGLMs( counts_data, count ~ exp_factor )
	fit1 = fitNbinomGLMs( counts_data, count ~ exp_factor + condition )
	glm_pvals = nbinomGLMTest( fit1, fit0 )

	results_table$pval <- glm_pvals
	results_table$log2FoldChange <- fit1[, grep("condition", colnames(fit1))]
	results_table$foldChange <- 2^fit1[, grep("condition", colnames(fit1))]
} 

# Adjust p-values
results_table$padj <- p.adjust(results_table$pval, method=p.value.adjustment.method)

# Merge with original data table and keep significant DEGs
significant_table <- cbind(dat, results_table)[results_table$padj <= p.value.cutoff, ]

# Remove rows where adjusted p-value is NA
significant_table <- significant_table[! (is.na(significant_table$padj)), ]

# Order by p-value
significant_table <- significant_table[ order(significant_table$padj), ] 

# Output significant DEGs
if (dim(significant_table)[1] > 0) {
	ndat <- ncol(dat)
	nmax <- ncol(significant_table)
	write.table(cbind(significant_table[,1:ndat], round(significant_table[, (ndat+1):(nmax-2)], digits=2), format(significant_table[, (nmax-1):nmax], digits=4, scientific=T)), file="de-list-deseq.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# Create a template output table for plotting. 
output_table <- NULL
colnames(results_table) <- gsub("\\..*$", "", colnames(results_table))
for(i in grep("baseMean$", colnames(results_table))) {
	col_size <- grep("padj", colnames(results_table))[1] - grep("baseMean", colnames(results_table))[1]
	output_table <- rbind(output_table, cbind(dat, results_table[, (i:(i+col_size))]))
}
rownames(output_table) <- make.names(rep(rownames(results_table), length(grep("baseMean$", colnames(results_table)))), unique=T)


# If genomic coordinates are present, output a sorted BED file for genome browser visualization and region matching tools
if("chr" %in% colnames(dat)) {
	if (dim(significant_table)[1] > 0) {
		bed_output <- significant_table [,c("chr","start","end")]
		bed_output <- cbind(bed_output,name=rownames(significant_table))
		bed_output <- cbind(bed_output, score=significant_table[,"log2FoldChange"])
		bed_output <- sort.bed(bed_output)
		write.table(bed_output, file="de-list-deseq.bed", sep="\t", row.names=F, col.names=F, quote=F)
	}
}

# Make dispersion plot
pdf(file="02-dispersion-plot-deseq.pdf")
plotDispEsts(counts_data, main="Dispersion plot", cex=0.2)
legend(x="topright", legend="fitted dispersion", col="red", cex=1, pch="-")
dev.off()

# Make histogram of p-values with overlaid significance cutoff and uniform distribution.
pdf (file="03-p-value-plot-deseq.pdf")
hist(output_table$pval, breaks=100, col="blue", border="slateblue", freq=FALSE, main="P-value distribution", xlab="p-value", ylab="proportion (%)")
hist(output_table$padj, breaks=100, col="red", border="slateblue", add=TRUE, freq=FALSE)
abline(h=1, lwd=2, lty=2, col="black")
abline(v=p.value.cutoff, lwd=2, lty=2, col="green")
legend (x="topright", legend=c("p-values","adjusted p-values", "uniform distribution", "significance cutoff"), col=c("blue","red","black","green"), cex=1, pch=15)
dev.off()

# Define function for making MA-plot of significant findings
plotDE <- function(res)
	plot(res$baseMean, res$log2FoldChange,
			log="x", pch=20, cex=.25, col = ifelse( res$padj <= p.value.cutoff, "red", "black"),
			main="MA plot", xlab="mean counts", ylab="log2(fold change)") 

# Make MA-plot
pdf(file="01-ma-plot-deseq.pdf")
plotDE(unique(output_table))
legend (x="topleft", legend=c("significant","not significant"), col=c("red","black"), cex=1, pch=19)
abline(h = c(-1, 0, 1), col = c("dodgerblue", "darkgreen", "dodgerblue"), lwd = 2)
dev.off()

# Join the PDFs 
system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=deseq_report.pdf *.pdf")

# EOF
