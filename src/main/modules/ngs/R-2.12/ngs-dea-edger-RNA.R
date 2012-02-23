# TOOL ngs-dea-edger-RNA.R: "Differential expression analysis using edgeR" (This tool will perform an analysis for differentially expressed sequences using the R implementation of the edge algorithm.)
# INPUT data.tsv TYPE GENERIC
# INPUT phenodata.tsv TYPE GENERIC
# OUTPUT OPTIONAL de-list-edger.tsv
# OUTPUT OPTIONAL de-list-edger.bed
# OUTPUT OPTIONAL ma-plot-raw-edger.pdf
# OUTPUT OPTIONAL ma-plot-normalized-edger.pdf
# OUTPUT OPTIONAL ma-plot-significant-edger.pdf
# OUTPUT OPTIONAL mds-plot-edger.pdf
# OUTPUT OPTIONAL edger-log.txt
# OUTPUT OPTIONAL p-value-plot-edger.pdf
# PARAMETER column: "Column describing groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER normalization: "Apply normalization" TYPE [yes, no] DEFAULT yes (If enabled, a normalization factor based on the trimmed mean of M-values \(TMM\) is performed to reduce the effect from sequencing biases.)
# PARAMETER dispersion_method: "Dispersion method" TYPE [common, tagwise] DEFAULT tagwise (The dispersion of counts for any given sequence can either be estimated based on the actual counts in the sample data set or be moderated across a selection of sequences with similar count numbers. The latter option, which is set by default, typically yields higher sensitivity and specificity. Note that when no biological replicates are available common dispersion is used regardless of the setting.)
# PARAMETER dispersion_estimate:"Dispersion estimate" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.1 (The value to use for estimating the common dispersion when no replicates are available.) 
# PARAMETER p_value_adjustment_method: "Multiple testing correction" TYPE [none, Bonferroni, Holm, Hochberg, BH, BY] DEFAULT BH (Multiple testing correction method.)
# PARAMETER p_value_threshold: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (The cutoff for statistical significance.)
# PARAMETER image_width: "Plot width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image_height: "Plot height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


############################################################
#                                                          #
# Analaysis workflow using edgeR for normalization and     #
# statistical testing for finding differentially expressed #
# sequence tags                                            #
#                                                          #
# MG, 11.6.2011                                            #
# updated, MG, 23.08.2011, to include library size from    #
# phenodata file                                           #
# updated MG, 30.01.2012 to allow analysis without         #
# biological replicates                                    #
# updated MG, 22.02.2012, prettified plots, added p-value  #
# distribution plot                                        #
#                                                          #
############################################################

# Loads the libraries
library(edgeR)

# Simplify variable names
w <- image_width
h <- image_height

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
if (number_samples == 2) dispersion_method <- "common" 

# Create a DGEList
# Notice that Library size is calculated from column totals if no library size
# exist in the phenodata file
if (estimate_lib_size) {
	dge_list <- DGEList (count=dat2, group=groups)
} else {
	dge_list <- DGEList (count=dat2, group=groups, lib.size=lib_size)
}

# Check lib size totals
# dge_list$samples

# Calculate normalization factors
if (normalization == "yes") {
	dge_list <- calcNormFactors(dge_list) 
}

# Produce MDS plot of normazied data
# NOTE: only possible when there are more than 2 samples in total
if (number_samples > 2) {
	pdf(file="mds-plot-edger.pdf", width=w/72, height=h/72)
	sample_colors <-  ifelse (dge_list$samples$group==group_levels[1], 1, 2)
	plotMDS.dge(dge_list, main="MDS Plot", col=sample_colors)
	legend(x="topleft", legend = group_levels,col=c(1,2), pch=19)
	dev.off()
}

# MA-plot comparison before and after normalization
pdf(file="ma-plot-raw-edger.pdf", width=w/72, height=h/72)
maPlot(dge_list$counts[, 1], dge_list$counts[, 2], normalize = FALSE, pch = 19,
		cex = 0.4, ylim = c(-8, 8))
grid(col = "blue")
title("Raw counts")
abline(h = log2(dge_list$samples$norm.factors[2]/dge_list$samples$norm.factors[1]),
		col = "red", lwd = 2)
abline(h = 0, col = "darkgreen", lwd = 1)
legend (x="topleft", legend=c("not epressed in one condition","expressed in both conditions"), col=c("orange","black"),
		cex=1, pch=19)
dev.off()

if (normalization == "yes") {
	pdf(file="ma-plot-normalized-edger.pdf", width=w/72, height=h/72)
	eff.libsize <- dge_list$samples$lib.size * dge_list$samples$norm.factors
	maPlot(dge_list$counts[, 1]/eff.libsize[1], dge_list$counts[, 2]/eff.libsize[2],
			normalize = TRUE, pch = 19, cex = 0.4, ylim = c(-8, 8))
	grid(col = "blue")
	title("Normalized counts")
	legend (x="topleft", legend=c("not epressed in one condition","expressed in both conditions"), col=c("orange","black"),
			cex=1, pch=19)
	abline(h = 0, col = "darkgreen", lwd = 1)
	dev.off()
}

###################################
# Anlysis using common dispersion #
###################################

if (dispersion_method == "common") {
	# Calculate common dispersion
	dge_list <- estimateCommonDisp(dge_list)
	
	# Statistical testing
	if (number_samples != 2) stat_test <- exactTest(dge_list, common.disp=TRUE) 
	if (number_samples == 2) stat_test <- exactTest(dge_list, common.disp=TRUE, dispersion=dispersion_estimate)
	
	
	# Extract results in a nice-looking table
	number_tags <- dim (dge_list$counts) [1]
	results_table <- topTags (stat_test, n=number_tags, adjust.method=p_value_adjustment_method, sort.by="p.value")
	results_table <- results_table$table
	
	# Extract the significant tags based on adjusted p-value cutoff
	cutoff <- p_value_threshold
	significant_results <- results_table[results_table$FDR<cutoff,]
	
	# Make an MA-plot displaying the significant reads
	pdf(file="ma-plot-significant-edger.pdf", width=w/72, height=h/72)	
	significant_indices <- rownames (significant_results)
	plotSmear(dge_list, de.tags = significant_indices, main = "MA plot for significantly\ndifferentially expressed sequence tags")
	abline(h = c(-1, 0, 1), col = c("dodgerblue", "darkgreen", "dodgerblue"), lwd = 2)
	legend (x="topleft", legend=c("significant","not significant"), col=c("red","black"),
			cex=1, pch=19)
	dev.off()
}

##############################################
# Anlysis using moderated tagwise dispersion #
##############################################

if (dispersion_method == "tagwise") {
	# Calculate the tagwise dispersion
	number_moderating_tags <- 10
	dge_list <- estimateCommonDisp(dge_list)
	dge_list <- estimateTagwiseDisp(dge_list, prior.n = number_moderating_tags)
	
	# Statistical testing
	stat_test <- exactTest(dge_list, common.disp = FALSE)
	
	# Extract results in a nice-looking table
	number_tags <- dim (dge_list$counts) [1]
	results_table <- topTags (stat_test, n=number_tags, adjust.method=p_value_adjustment_method, sort.by="p.value")
	results_table <- results_table$table
	
	# Extract the significant tags based on adjusted p-value cutoff
	cutoff <- p_value_threshold
	significant_results <- results_table[results_table$FDR<cutoff,]
	
	# Make an MA-plot displaying the significant reads
	pdf(file="ma-plot-significant-edger.pdf", width=w/72, height=h/72)	
	significant_indices <- rownames (significant_results)
	plotSmear(dge_list, de.tags = significant_indices, main = "MA plot for significantly\ndifferentially expressed sequence tags")
	abline(h = c(-1, 0, 1), col = c("dodgerblue", "darkgreen", "dodgerblue"), lwd = 2)
	legend (x="topleft", legend=c("significant features","not significant"), col=c("red","black"),
			cex=1, pch=19)
	dev.off()
}

# Create a tbale with the original counts per sample together with the statistical tests results
# ready for output in Chipster
# If there are no significant results return a message
if (dim(significant_results)[1] > 0) {
	output_table <- data.frame (dat[significant_indices,], significant_results)
}

# Output the table
if (dim(significant_results)[1] > 0) {
	write.table(output_table, file="de-list-edger.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# Also output a bed graph file for visualization and region matching tools
if (dim(significant_results)[1] > 0) {
	empty_column <- character(length(significant_indices))
	bed_output <- output_table [,c("chr","start","end")]
	bed_output <- cbind(bed_output,empty_column)
	bed_output <- cbind(bed_output, output_table[,"logFC"])
	write.table(bed_output, file="de-list-edger.bed", sep="\t", row.names=F, col.names=F, quote=F)
}

# Output a message if no significant genes are found
if (dim(significant_results)[1] == 0) {
	cat("No statistically significantly expressed sequences were found. Try again with a less stringent p-value cut-off or multiple testing correction method.", file="edger-log.txt")
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
