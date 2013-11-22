# TOOL ngs-dea-edger-RNA.R: "Differential expression using edgeR" (Differential gene expression analysis using the exact test of the edgeR Bioconductor package. You can create the input count table and phenodata file using the tool \"Utilities - Define NGS experiment\". You should set the filtering parameter to the number of samples in your smallest experimental group. Please note that this tool is suitable only for two group comparisons. For multifactor experiments please use the tool \"Differential expression using edgeR for multivariate experiments\".)
# INPUT data.tsv TYPE GENERIC
# INPUT phenodata.tsv TYPE GENERIC
# OUTPUT OPTIONAL de-list-edger.tsv
# OUTPUT OPTIONAL de-list-edger.bed
# OUTPUT OPTIONAL ma-plot-edger.pdf
# OUTPUT OPTIONAL mds-plot-edger.pdf
# OUTPUT OPTIONAL edger-log.txt
# OUTPUT OPTIONAL p-value-plot-edger.pdf
# OUTPUT OPTIONAL dispersion-edger.pdf
# PARAMETER column: "Column describing groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test.)
# PARAMETER filter: "Filter out genes which don't have counts in at least this many samples" TYPE INTEGER FROM 1 TO 1000 DEFAULT 2 (Analyze only genes which have at least 5 counts in at least this many samples. You should set this to the number of samples in your smallest experimental group.)
# PARAMETER OPTIONAL p_value_threshold: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (The cutoff for adjusted p-values.)
# PARAMETER OPTIONAL p_value_adjustment_method: "Multiple testing correction" TYPE [none, Bonferroni, Holm, Hochberg, BH, BY] DEFAULT BH (Multiple testing correction method.)
# PARAMETER OPTIONAL dispersion_method: "Dispersion method" TYPE [common, tagwise] DEFAULT tagwise (The dispersion of counts for a gene can be moderated across several genes with similar count numbers. This default tagwise option typically yields higher sensitivity and specificity. The option Common estimates one value which is then used for all the genes. Common dispersion is used regardless of the setting if no biological replicates are available.)
# PARAMETER OPTIONAL dispersion_estimate:"Dispersion value used if no replicates are available" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.1 (The value to use for the common dispersion when no replicates are available.) 
# PARAMETER OPTIONAL normalization: "Apply TMM normalization" TYPE [yes, no] DEFAULT yes (Should normalization based on the trimmed mean of M-values \(TMM\) be performed to reduce the RNA composition effect.)
# PARAMETER OPTIONAL w: "Plot width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted image)
# PARAMETER OPTIONAL h: "Plot height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted image)


# MG 11.6.2011                                            
# MG 23.8.2011, included library size from phenodata file                                           
# MG 30.1.2012, allowed analysis without biological replicates                                    
# MG 22.2.2012, prettified plots, added p-value  distribution plot
# EK 6.5.2012, clarified wording
# AMS 5.10.2012, added sorting to BED
# EK 30.4.2013, changes to descriptions, made genomic location info optional so that external count tables can be used
# EK 2.5.2013, added dispersion plot and filtering based on counts, disabled extra MA plots
# EK 5.5.2013, modified filtering based on counts, removed fixed prior.n
# EK 19.11.2013, updated to edgeR 3.4.0 and enabled trended dispersion. Filtering det by default to 1.

# OUTPUT OPTIONAL ma-plot-raw-edger.pdf
# OUTPUT OPTIONAL ma-plot-normalized-edger.pdf

# Loads the libraries
library(edgeR)

# Loads the count data
file <- c("data.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)

# Extracts expression value columns
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
if (number_samples == 2) dispersion_method <- "common" 

# Create a DGEList. If no library size exist in the phenodata file, it is calculated from column totals 
if (estimate_lib_size) {
	dge_list <- DGEList (count=dat2, group=groups)
} else {
	dge_list <- DGEList (count=dat2, group=groups, lib.size=lib_size)
}

# filter out genes which have less than 5 counts in user-defined number of samples 
if (filter > 0) {
	keep <- rowSums(dge_list$counts>5) >= filter
	dge_list <- dge_list[keep,]
	dge_list$lib.size <- colSums(dge_list$counts)
}

# to be considered later: filter based on counts per million
# m <- sweep(dge_list$counts, 2, 1e6 / dge_list$samples$lib.size, `*`)
# ridx <- rowSums(m > 1) >= 3
# dge_list <- dge_list[ridx,]

# Calculate normalization factors
if (normalization == "yes") {
	dge_list <- calcNormFactors(dge_list) 
}


# Produce MDS plot of normalized data (possible only when there are more than 2 samples)
if (number_samples > 2) {
	pdf(file="mds-plot-edger.pdf", width=w/72, height=h/72)
	sample_colors <-  ifelse (dge_list$samples$group==group_levels[1], 1, 2)
	plotMDS(dge_list, main="MDS Plot", col=sample_colors)
	legend(x="topleft", legend = group_levels,col=c(1,2), pch=19)
	dev.off()
}

# MA plots before and after normalization (currently disabled output)
pdf(file="ma-plot-raw-edger.pdf", width=w/72, height=h/72)
maPlot(dge_list$counts[, 1], dge_list$counts[, 2], normalize = FALSE, pch = 19, cex = 0.4, ylim = c(-8, 8))
grid(col = "blue")
title("Raw counts")
abline(h = log2(dge_list$samples$norm.factors[2]/dge_list$samples$norm.factors[1]), col = "red", lwd = 2)
abline(h = 0, col = "darkgreen", lwd = 1)
legend (x="topleft", legend=c("not expressed in one condition","expressed in both conditions"), col=c("orange","black"), cex=1, pch=19)
dev.off()

if (normalization == "yes") {
	pdf(file="ma-plot-normalized-edger.pdf", width=w/72, height=h/72)
	eff.libsize <- dge_list$samples$lib.size * dge_list$samples$norm.factors
	maPlot(dge_list$counts[, 1]/eff.libsize[1], dge_list$counts[, 2]/eff.libsize[2],
			normalize = TRUE, pch = 19, cex = 0.4, ylim = c(-8, 8))
	grid(col = "blue")
	title("Normalized counts")
	legend (x="topleft", legend=c("not expressed in one condition","expressed in both conditions"), col=c("orange","black"), cex=1, pch=19)
	abline(h = 0, col = "darkgreen", lwd = 1)
	dev.off()
}

# Analysis using common dispersion
if (dispersion_method == "common") {
	# Calculate common dispersion
	dge_list <- estimateCommonDisp(dge_list)
	
	# Statistical testing
	if (number_samples != 2) stat_test <- exactTest(dge_list) 
	if (number_samples == 2) stat_test <- exactTest(dge_list, dispersion=dispersion_estimate)
}	

# Analysis using moderated tagwise dispersion
if (dispersion_method == "tagwise") {
	# Calculate the tagwise dispersion
	dge_list <- estimateCommonDisp(dge_list)
	dge_list <- estimateTrendedDisp(dge_list)
	dge_list <- estimateTagwiseDisp(dge_list)
	# Statistical testing
	stat_test <- exactTest(dge_list)
	
# Dispersion plot
	pdf(file="dispersion-edger.pdf", width=w/72, height=h/72)
	plotBCV(dge_list, main="Biological coefficient of variation")
	dev.off()
}

# Extract results in a nice-looking table
number_tags <- dim (dge_list$counts) [1]
results_table <- topTags (stat_test, n=number_tags, adjust.method=p_value_adjustment_method, sort.by="p.value")
results_table <- results_table$table

# Extract the significant tags based on adjusted p-value cutoff
significant_results <- results_table[results_table$FDR<p_value_threshold,]

# Make an MA-plot displaying the significant reads
pdf(file="ma-plot-edger.pdf", width=w/72, height=h/72)	
significant_indices <- rownames (significant_results)
plotSmear(dge_list, de.tags = significant_indices, main = "MA plot")
abline(h = c(-1, 0, 1), col = c("dodgerblue", "darkgreen", "dodgerblue"), lwd = 2)
legend (x="topleft", legend=c("significantly differentially expressed","not significant"), col=c("red","black"),cex=1, pch=19)
dev.off()

# If significant results are found, create an output table with the original counts per sample together with the statistical tests results
# If genomic coordinates are present, output a sorted BED file for genome browser visualization and region matching tools
if (dim(significant_results)[1] > 0) {
	output_table <- data.frame (dat[significant_indices,], significant_results)
	write.table(output_table, file="de-list-edger.tsv", sep="\t", row.names=T, col.names=T, quote=F)
	these.colnames <- colnames(dat)
	if("chr" %in% these.colnames) {
		empty_column <- character(length(significant_indices))
		bed_output <- output_table [,c("chr","start","end")]
		bed_output <- cbind(bed_output,empty_column)
		bed_output <- cbind(bed_output, output_table[,"logFC"])
		source(file.path(chipster.common.path, "bed-utils.R"))
		bed_output <- sort.bed(bed_output)
		write.table(bed_output, file="de-list-edger.bed", sep="\t", row.names=F, col.names=F, quote=F)
	}	
}


# Output a message if no significant genes are found
if (dim(significant_results)[1] == 0) {
	cat("No statistically significantly expressed genes were found.", file="edger-log.txt")
}


# Make histogram of p-values with overlaid significance cutoff and uniform distribution
pdf (file="p-value-plot-edger.pdf")
hist(results_table$PValue, breaks=100, col="blue",
		border="slateblue", freq=FALSE,
		main="P-value distribution", xlab="p-value", ylab="proportion (%)")
hist(results_table$FDR, breaks=100, col="red",
		border="slateblue", add=TRUE, freq=FALSE)
abline(h=1, lwd=2, lty=2, col="black")
abline(v=p_value_threshold, lwd=2, lty=2, col="green")
legend (x="topright", legend=c("p-values","adjusted p-values", "uniform distribution", "significance cutoff"), col=c("blue","red","black","green"),
		cex=1, pch=15)
dev.off()

# EOF
