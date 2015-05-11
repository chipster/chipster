# TOOL deseq2.R: "Differential expression using DESeq2" (Differential expression analysis using the DESeq2 Bioconductor package. You can create the input count table and phenodata file using the tool \"Utilities - Define NGS experiment\". If using DESeq2 for finding differentially expressed genes between more than two experimental groups, note that output figures sum up information from all pairwise comparisons.)
# INPUT data.tsv TYPE GENERIC
# INPUT META phenodata.tsv TYPE GENERIC
# OUTPUT OPTIONAL de-list-deseq2.tsv
# OUTPUT OPTIONAL de-list-deseq2.bed
# OUTPUT OPTIONAL ma-plot-deseq2.pdf
# OUTPUT OPTIONAL dispersion-plot-deseq2.pdf
# OUTPUT OPTIONAL p-value-plot-deseq2.pdf
# PARAMETER column: "Column describing groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test.)
# PARAMETER ad_factor: "Column describing additional experimental factor" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing an additional experimental factor. If given, p-values in the output table are from a likelihood ratio test of a model including the experimental groups and experimental factor vs a model which only includes the experimental factor.)
# PARAMETER OPTIONAL normalization: "Apply normalization" TYPE [yes, no] DEFAULT yes (Should effective library size be estimated. This corrects for RNA composition bias. Note that if you have supplied library size in phenodata, size factors are calculated based on the library size total, and composition bias is not corrected.)
# PARAMETER OPTIONAL dispersion_estimate:"Dispersion estimation method" TYPE [parametric: "parametric", local: "local"] DEFAULT parametric (Dispersion can be estimated using a local fit or a two-coefficient parametric model. You should use local fit if there are no biological replicates.)
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

# Loads the correct library
source(file.path(chipster.common.path, "bed-utils.R"))
library(DESeq2)

# Loads the counts data and extract expression values
file <- c("data.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)
annotations <- dat[,-grep("chip", names(dat))]
dat2 <- dat[,grep("chip", names(dat))]

# Test needs a parameter "groups" that specifies experimental conditions
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
groups <- as.character (phenodata[,pmatch(column,colnames(phenodata))])
group_levels <- levels(as.factor(groups))

# Read addittional factors for GLM test and DESeq2. Construct design matrix from this 
exp_factor <- NULL
if(ad_factor != "EMPTY") {
	exp_factor <- as.character (phenodata[,pmatch(ad_factor,colnames(phenodata))])
	design <- data.frame(condition=as.factor(groups), exp_factor=exp_factor)
	rownames(design) <- colnames(dat2)
}

# Create a counts data object
if (ad_factor == "EMPTY") {
	counts_data <- DESeqDataSetFromMatrix(countData=dat2, colData=data.frame(condition=groups), design = ~ condition)
} else if (ad_factor != "EMPTY") {
	counts_data <- DESeqDataSetFromMatrix(countData=dat2, colData=design, design = ~ exp_factor + condition)
}

# Calculate size factors based on estimated library size, unless it is given in phenodata
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
counts_data <- estimateDispersions(counts_data, fitType=dispersion_estimate)	


#vector / variable that holds comparison names
results_name <- NULL

# Calculate statistic for differential expression
if (length(unique(groups)) > 2) {
	# If using DESeq >= 1.3, no need to set betaPrior=FALSE
	test_results <- nbinomWaldTest(counts_data, betaPrior = FALSE)		
	results_table <- NULL
	for (i in levels(colData(test_results)$condition)[-(length(levels(colData(test_results)$condition)))]) {
		for (j in levels(colData(test_results)$condition)[-(1:i)]) {
			pairwise_results <- as.data.frame(results(test_results, contrast=c("condition",i,j)))
			pairwise_results$padj <- p.adjust(pairwise_results$pval, method=p.value.adjustment.method)
			if(is.null(results_table)) results_table <- pairwise_results else  results_table <- cbind(results_table, pairwise_results)
			results_name <- c(results_name, paste(i,"_vs_", j, sep=""))
		}
	}
	colnames(results_table) <- paste(colnames(results_table), rep(results_name,each=6), sep=".")
} else if (length(unique(groups)) == 2 && length(unique(exp_factor)) > 2) {
	# If using DESeq >= 1.3, no need to set betaPrior=FALSE
	results_table <- results(nbinomWaldTest(counts_data, betaPrior = FALSE))
} else {
	results_table <- results(nbinomWaldTest(counts_data))
}

# Merge with original data table and keep significant DEGs
if(length(unique(groups)) == 2) {
	results_table$padj <- p.adjust(results_table$pval, method=p.value.adjustment.method)
	significant_table <- cbind(dat, results_table)[results_table$padj <= p.value.cutoff, ]
	significant_table <- significant_table[! (is.na(significant_table$padj)), ]
	significant_table <- significant_table[ order(significant_table$padj), ] 
} else {
	min_padj <- apply(results_table[, grep("padj", colnames(results_table))], 1, min)
	significant_table <- cbind(dat, results_table, min_padj=min_padj)
	significant_table <- significant_table[ (significant_table$min_padj <=  p.value.cutoff), ]
	significant_table <- significant_table[! (is.na(significant_table$min_padj)), ]
	significant_table <- significant_table[ order(significant_table$min_padj), ] 
	significant_table <- significant_table[, -grep("min_padj", colnames(significant_table))]
}

# Output significant DEGs
if (dim(significant_table)[1] > 0) {
	ndat <- ncol(dat)
	nmax <- ncol(significant_table)
	write.table(cbind(significant_table[,1:ndat], round(significant_table[, (ndat+1):(nmax-2)], digits=2), format(significant_table[, (nmax-1):nmax], digits=4, scientific=T)), file="de-list-deseq2.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

#Create a template output table for plotting. If having N genes and 3 comparisons, this conversion results in a data matrix that has Nx3 rows 
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
		bed_output <- output_table[,c("chr","start","end")]
		if(is.null(results_name)) {
			#peak_names <- character(nrow(output_table))
			peak_names <- rownames(results_table)
		} else {
			#peak_names <- paste(rep(results_name, each=nrow(results_table)), "_peak", 1:nrow(results_table), sep="")
			peak_names <- paste(rep(results_name, each=nrow(results_table)), rownames(results_table), sep="")	
		}
		bed_output <- cbind(bed_output, name=peak_names)							#name
		bed_output <- cbind(bed_output, score=output_table[, "log2FoldChange"])		#score
		bed_output <- bed_output[(output_table$padj <= p.value.cutoff & (! (is.na(output_table$padj)))), ]

		bed_output <- sort.bed(bed_output)
		write.table(bed_output, file="de-list-deseq2.bed", sep="\t", row.names=F, col.names=F, quote=F)
	}
}

# Make dispersion plot
pdf(file="dispersion-plot-deseq2.pdf")
plotDispEsts(counts_data, main="Dispersion plot", cex=0.2)
legend(x="topright", legend="fitted dispersion", col="red", cex=1, pch="-")
dev.off()

# Make histogram of p-values with overlaid significance cutoff and uniform distribution. When more than two groups, min.pvalue is taken over all comparisons for genes
pdf (file="p-value-plot-deseq2.pdf")
hist(output_table$pval, breaks=100, col="blue", border="slateblue", freq=FALSE, main="P-value distribution", xlab="p-value", ylab="proportion (%)")
hist(output_table$padj, breaks=100, col="red", border="slateblue", add=TRUE, freq=FALSE)
abline(h=1, lwd=2, lty=2, col="black")
abline(v=p.value.cutoff, lwd=2, lty=2, col="green")
legend (x="topright", legend=c("p-values","adjusted p-values", "uniform distribution", "significance cutoff"), col=c("blue","red","black","green"), cex=1, pch=15)
dev.off()

# Define function for making MA-plot of significant findings
plotDE <- function(res)
	plot(res$baseMean, res$log2FoldChange,
			log="x", pch=20, cex=.25, col = ifelse( res$padj < p.value.cutoff, "red", "black"),
			main="MA plot", xlab="mean counts", ylab="log2(fold change)") 

# Make MA-plot
pdf(file="ma-plot-deseq2.pdf")
plotDE(unique(output_table))
legend (x="topleft", legend=c("significant","not significant"), col=c("red","black"), cex=1, pch=19)
abline(h = c(-1, 0, 1), col = c("dodgerblue", "darkgreen", "dodgerblue"), lwd = 2)
dev.off()

# EOF
