# TOOL deseq2.R: "Differential expression using DESeq2" (Differential expression analysis using the DESeq2 Bioconductor package. You can create the input count table and phenodata file using the tool \"Utilities - Define NGS experiment\". If you have more than two experimental groups, note that the output figures sum up information from all pairwise comparisons.)
# INPUT data.tsv: "Count table" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata file" TYPE GENERIC
# OUTPUT OPTIONAL de-list-deseq2.tsv
# OUTPUT OPTIONAL summary.txt
# OUTPUT OPTIONAL deseq2_report.pdf
# PARAMETER column: "Column describing groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test.)
# PARAMETER OPTIONAL ad_factor: "Column describing additional experimental factor" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing an additional experimental factor. If given, p-values in the output table are from a likelihood ratio test of a model including the experimental groups and experimental factor, vs a model which only includes the experimental factor.)
# PARAMETER OPTIONAL p.value.cutoff: "Cutoff for the adjusted P-value" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (The cutoff for Benjamini-Hochberg adjusted p-value.)


# MK 15.04.2014, added the possibility to use DESeq2 in dea-deseq.R 
# AMS 17.06.2014, split the DESeq2 part to a separate tool
# EK 1.7.2014, clarified the script before moving it to production, and fixed a bug that disabled DESeq2's automatic independent filtering 
# EK 9.2.2015, updated to R3.1.2, changed the MA plot, added summary
# AMS 7.4.2015, Join pdf outputs to one

#column <-"group"
#ad_factor<-"EMPTY"
#p.value.cutoff<-0.05

# Load the library
source(file.path(chipster.common.path, "bed-utils.R"))
library(DESeq2)

# Load the counts data and extract expression value columns
dat <- read.table("data.tsv", header=T, sep="\t", row.names=1)
dat2 <- dat[,grep("chip", names(dat))]

# Get the experimental group information from phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
groups <- as.character (phenodata[,pmatch(column,colnames(phenodata))])
group_levels <- levels(as.factor(groups))

# Read the additional experimental factor from phenodata and construct a design matrix from this information
exp_factor <- NULL
if(ad_factor != "EMPTY") {
	exp_factor <- as.character (phenodata[,pmatch(ad_factor,colnames(phenodata))])
	design <- data.frame(condition=as.factor(groups), exp_factor=exp_factor)
	rownames(design) <- colnames(dat2)
}

# Create a DESeqDataSet object
if (ad_factor == "EMPTY") {
	dds <- DESeqDataSetFromMatrix(countData=dat2, colData=data.frame(condition=groups), design = ~ condition)
} else if (ad_factor != "EMPTY") {
	dds <- DESeqDataSetFromMatrix(countData=dat2, colData=design, design = ~ exp_factor + condition)
}
	
# Vector / variable that holds comparison names
results_name <- NULL 

# Calculate statistic for differential expression, merge with original data table, keep significant DEGs, remove NAs and sort by FDR. If there are more than 2 groups, get pairwise results for each comparison.
if (length(unique(groups)) == 2) {
	dds <- DESeq(dds)
	res <- results(nbinomWaldTest(dds))
	sig <- cbind(dat, res)[res$padj <= p.value.cutoff, ]
	sig <- sig[! (is.na(sig$padj)), ]
	sig <- sig[ order(sig$padj), ]
	# Open pdf file for output
	pdf(file="deseq2_report.pdf") 
	plotMA(dds,alpha=p.value.cutoff,main=c("DESeq2 MA-plot, FDR =", p.value.cutoff),ylim=c(-2,2))
	sink("summary.txt")
	summary(res, alpha=p.value.cutoff)
	sink()
	
} else if (length(unique(groups)) > 2){
	dds <- estimateSizeFactors(dds)
	dds <- estimateDispersions(dds)
	test_results <- nbinomWaldTest(dds)
	res <- NULL
	for (i in levels(colData(test_results)$condition)[-(length(levels(colData(test_results)$condition)))]) {
		for (j in levels(colData(test_results)$condition)[-(1:i)]) {
			pairwise_results <- as.data.frame(results(test_results, contrast=c("condition",i,j)))
			if(is.null(res)) res <- pairwise_results else  res <- cbind(res, pairwise_results)
			results_name <- c(results_name, paste(i,"_vs_", j, sep=""))
		}
	}
	colnames(res) <- paste(colnames(res), rep(results_name,each=6), sep=".")
	min_padj <- apply(res[, grep("padj", colnames(res))], 1, min)
	sig <- cbind(dat, res, min_padj=min_padj)
	sig <- sig[ (sig$min_padj <=  p.value.cutoff), ]
	sig <- sig[! (is.na(sig$min_padj)), ]
	sig <- sig[ order(sig$min_padj), ] 
	sig <- sig[, -grep("min_padj", colnames(sig))]
} 


# Output significant DEGs
if (dim(sig)[1] > 0) {
	ndat <- ncol(dat)
	nmax <- ncol(sig)
	write.table(cbind(sig[,1:ndat], round(sig[, (ndat+1):(nmax-2)], digits=2), format(sig[, (nmax-1):nmax], digits=4, scientific=T)), file="de-list-deseq2.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# Create a template output table for plotting. If having N genes and 3 comparisons, this conversion results in a data matrix that has Nx3 rows 
output_table <- NULL
colnames(res) <- gsub("\\..*$", "", colnames(res))
for(i in grep("baseMean$", colnames(res))) {
	col_size <- grep("padj", colnames(res))[1] - grep("baseMean", colnames(res))[1]
	output_table <- rbind(output_table, cbind(dat, res[, (i:(i+col_size))]))
}
rownames(output_table) <- make.names(rep(rownames(res), length(grep("baseMean$", colnames(res)))), unique=T)

# If genomic coordinates are present, output a sorted BED file for genome browser visualization and region matching tools
if("chr" %in% colnames(dat)) {
	if (dim(sig)[1] > 0) {
		bed <- output_table[,c("chr","start","end")]
		if(is.null(results_name)) {
			gene_names <- rownames(res)
		} else {
			gene_names <- paste(rep(results_name, each=nrow(res)), rownames(res), sep="")	
		}
		bed <- cbind(bed, name=gene_names)							#name
		bed <- cbind(bed, score=output_table[, "log2FoldChange"])		#score
		bed <- bed[(output_table$padj <= p.value.cutoff & (! (is.na(output_table$padj)))), ]
		bed <- sort.bed(bed)
		write.table(bed, file="de-list-deseq2.bed", sep="\t", row.names=F, col.names=F, quote=F)
	}
}

# Make dispersion plot
plotDispEsts(dds, main="Dispersion plot", cex=0.2)
legend(x="topright", legend="fitted dispersion", col="red", cex=1, pch="-")

# Make histogram of p-values with overlaid significance cutoff. When more than two groups, min.pvalue is taken over all comparisons for genes
hist(output_table$pval, breaks=100, col="blue", border="slateblue", freq=FALSE, main="P-value distribution", xlab="p-value", ylab="proportion (%)")
hist(output_table$padj, breaks=100, col="red", border="slateblue", add=TRUE, freq=FALSE)
abline(v=p.value.cutoff, lwd=2, lty=2, col="black")
legend (x="topright", legend=c("p-values","adjusted p-values", "significance cutoff"), col=c("blue","red","black"), cex=1, pch=15)
# Close pdf
dev.off()

# MA-plot when there are more than 2 groups. Define function for making MA-plot.
# plotDE <- function(res)
#	plot(res$baseMean, res$log2FoldChange,
#			log="x", pch=20, cex=.25, col = ifelse( res$padj < p.value.cutoff, "red", "black"),
#			main="MA plot", xlab="mean counts", ylab="log2(fold change)") 
# Make MA-plot
# pdf(file="ma-plot-deseq2.pdf")
# plotDE(unique(output_table))
# legend (x="topleft", legend=c("significant","not significant"), col=c("red","black"), cex=1, pch=19)
# abline(h = c(-1, 0, 1), col = c("dodgerblue", "darkgreen", "dodgerblue"), lwd = 2)
# dev.off()

# EOF
