# TOOL deseq2.R: "Differential expression using DESeq2" (Differential expression analysis using the DESeq2 Bioconductor package. You can create the input count table and phenodata file using the tool \"Utilities - Define NGS experiment\". If you have more than two experimental groups, note that the output figures sum up information from all pairwise comparisons.)
# INPUT data.tsv TYPE GENERIC
# INPUT phenodata.tsv TYPE GENERIC
# OUTPUT OPTIONAL de-list-deseq2.tsv
# OUTPUT OPTIONAL de-list-deseq2.bed
# OUTPUT OPTIONAL ma-plot-deseq2.pdf
# OUTPUT OPTIONAL dispersion-plot-deseq2.pdf
# OUTPUT OPTIONAL p-value-plot-deseq2.pdf
# PARAMETER column: "Column describing groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test.)
# PARAMETER OPTIONAL ad_factor: "Column describing additional experimental factor" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing an additional experimental factor. If given, p-values in the output table are from a likelihood ratio test of a model including the experimental groups and experimental factor, vs a model which only includes the experimental factor.)
# PARAMETER OPTIONAL dispersion_estimate:"Dispersion estimation method" TYPE [parametric: "parametric", local: "local"] DEFAULT parametric (Dispersion can be estimated using a local fit or a two-coefficient gamma-family GLM. You should use local fit if there are no biological replicates.)
# PARAMETER OPTIONAL p.value.cutoff: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (The cutoff for adjusted p-value.)


# MK 15.04.2014, added the possibility to use DESeq2 in dea-deseq.R 
# AMS 17.06.2014, split the DESeq2 part to a separate tool
# EK 1.7.2014, clarified the script before moving the it to production and fixed the bug that disabled DESeq2's automatic independent filtering 


#column <-"group"
#ad_factor<-"EMPTY"
#dispersion_estimate<-"parametric"
#p.value.cutoff<-0.05

# Loads the correct library
source(file.path(chipster.common.path, "bed-utils.R"))
library(DESeq2)

# Loads the counts data and extract expression values
file <- c("data.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)
dat2 <- dat[,grep("chip", names(dat))]

# Get the experimental group information from the relevant phenodata column
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
groups <- as.character (phenodata[,pmatch(column,colnames(phenodata))])
group_levels <- levels(as.factor(groups))

# Read additional experimental factors from phenodata and construct design matrix from this information
exp_factor <- NULL
if(ad_factor != "EMPTY") {
	exp_factor <- as.character (phenodata[,pmatch(ad_factor,colnames(phenodata))])
	design <- data.frame(condition=as.factor(groups), exp_factor=exp_factor)
	rownames(design) <- colnames(dat2)
}

# Create a DESeqDataSet object from the counts file
if (ad_factor == "EMPTY") {
	counts_data <- DESeqDataSetFromMatrix(countData=dat2, colData=data.frame(condition=groups), design = ~ condition)
} else if (ad_factor != "EMPTY") {
	counts_data <- DESeqDataSetFromMatrix(countData=dat2, colData=design, design = ~ exp_factor + condition)
}

# Calculate size factors based on estimated library size
counts_data <- estimateSizeFactors(counts_data)

# Estimate dispersion values for each gene. Use fitType to control for parametric or local fit. 
counts_data <- estimateDispersions(counts_data, fitType=dispersion_estimate)	

# Vector / variable that holds comparison names
results_name <- NULL 

# Calculate statistic for differential expression, merge with original data table, keep significant DEGs, remove NAs and sort by FDR. If there are more than 2 groups, get pairwise results for each comparison.
if (length(unique(groups)) == 2) {
	results_table <- results(nbinomWaldTest(counts_data))
	significant_table <- cbind(dat, results_table)[results_table$padj <= p.value.cutoff, ]
	significant_table <- significant_table[! (is.na(significant_table$padj)), ]
	significant_table <- significant_table[ order(significant_table$padj), ]
} else if (length(unique(groups)) > 2){
	test_results <- nbinomWaldTest(counts_data)
	results_table <- NULL
	for (i in levels(colData(test_results)$condition)[-(length(levels(colData(test_results)$condition)))]) {
		for (j in levels(colData(test_results)$condition)[-(1:i)]) {
			pairwise_results <- as.data.frame(results(test_results, contrast=c("condition",i,j)))
			if(is.null(results_table)) results_table <- pairwise_results else  results_table <- cbind(results_table, pairwise_results)
			results_name <- c(results_name, paste(i,"_vs_", j, sep=""))
		}
	}
	colnames(results_table) <- paste(colnames(results_table), rep(results_name,each=6), sep=".")
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
			gene_names <- rownames(results_table)
		} else {
			gene_names <- paste(rep(results_name, each=nrow(results_table)), rownames(results_table), sep="")	
		}
		bed_output <- cbind(bed_output, name=gene_names)							#name
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

# Make histogram of p-values with overlaid significance cutoff. When more than two groups, min.pvalue is taken over all comparisons for genes
pdf (file="p-value-plot-deseq2.pdf")
hist(output_table$pval, breaks=100, col="blue", border="slateblue", freq=FALSE, main="P-value distribution", xlab="p-value", ylab="proportion (%)")
hist(output_table$padj, breaks=100, col="red", border="slateblue", add=TRUE, freq=FALSE)
abline(v=p.value.cutoff, lwd=2, lty=2, col="black")
legend (x="topright", legend=c("p-values","adjusted p-values", "significance cutoff"), col=c("blue","red","black"), cex=1, pch=15)
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
