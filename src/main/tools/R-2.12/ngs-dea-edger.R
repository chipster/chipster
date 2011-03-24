# TOOL "Statistics" / ngs-dea-edger.R: "Differential expression analysis using edgeR" (This tool will perform an analysis for differentially expressed sequences using the R implementation of the edge algorithm.)
# INPUT data.tsv TYPE GENERIC
# INPUT phenodata.tsv TYPE GENERIC
# OUTPUT de-list.tsv
# OUTPUT OPTIONAL ma-plot-raw-counts.png
# OUTPUT OPTIONAL ma-plot-normalized-counts.png
# OUTPUT OPTIONAL ma-plot-significant-counts.png
# OUTPUT OPTIONAL mds-plot.png
# PARAMETER column: "Column describing groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER normalization: "Apply normalization" TYPE [yes, no] DEFAULT yes (If enabled, a normalization factor based on the trimmed mean of M-values \(TMM\) is performed to reduce the effect from sequencing biases.)
# PARAMETER dispersion_estimate: "Dispersion estimate" TYPE [common, tagwise] DEFAULT tagwise (The dispersion of counts for any given sequence can either be estimated based on the actual counts in the sample data set or be moderated across a selection of sequences with similar count numbers. The latter option, which is set by default, typically yields higher sensitivity and specificity.)
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
# MG, 11.3.2011                                            #
# development version                                      #
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
annotations <- dat[,-grep("Sample", names(dat))]
dat2 <- dat[,grep("Sample", names(dat))]

# Test needs a parameter "groups" that specifies the grouping of the samples
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
groups <- as.character (phenodata[,pmatch(column,colnames(phenodata))])
group_levels <- levels(as.factor(groups))

# Sanity checks
if(length(unique(groups))==1 | length(unique(groups))>=3) {
	stop("You need to have exactly two groups to run this analysis")
}

# Create a DGEList
# Notice that Library size is calculated from column totals
dge_list <- DGEList (count=dat2, group=groups)

# Check lib size totals
# dge_list$samples

# Calculate normalization factors
if (normalization == "yes") {
	dge_list <- calcNormFactors(dge_list) 
}

# Produce MDS plot of normazied data
bitmap(file="mds-plot.png", width=w/72, height=h/72)
sample_colors <-  ifelse (dge_list$samples$group==group_levels[1], 1, 2)
plotMDS.dge(dge_list, main="MDS Plot", xlim=c(-2,1), col=sample_colors)
legend(x="topleft", legend = group_levels,col=c(1,2), pch=19)
dev.off()

# MA-plot comparison before and after normalization
bitmap(file="ma-plot-raw-counts.png", width=w/72, height=h/72)
maPlot(dge_list$counts[, 1], dge_list$counts[, 2], normalize = TRUE, pch = 19,
		cex = 0.4, ylim = c(-8, 8))
grid(col = "blue")
title("Raw counts")
abline(h = log2(dge_list$samples$norm.factors[2]/dge_list$samples$norm.factors[1]),
		col = "red", lwd = 4)
dev.off()

if (normalization == "yes") {
	bitmap(file="ma-plot-normalized-counts.png", width=w/72, height=h/72)
	eff.libsize <- dge_list$samples$lib.size * dge_list$samples$norm.factors
	maPlot(dge_list$counts[, 1]/eff.libsize[1], dge_list$counts[, 2]/eff.libsize[2],
			normalize = FALSE, pch = 19, cex = 0.4, ylim = c(-8, 8))
	grid(col = "blue")
	title("Normalized counts")
	dev.off()
}

###################################
# Anlysis using common dispersion #
###################################

if (dispersion_estimate == "common") {
	# Calculate common dispersion
	dge_list <- estimateCommonDisp(dge_list)
	
	# Statistical testing
	stat_test <- exactTest(dge_list)
	
	# Extract results in a nice-looking table
	number_tags <- dim (dge_list$counts) [1]
	results_table <- topTags (stat_test, n=number_tags, adjust.method=p_value_adjustment_method, sort.by="p.value")
	results_table <- results_table$table
	
	# Extract the significant tags based on adjusted p-value cutoff
	cutoff <- p_value_threshold
	significant_results <- results_table[results_table$FDR<cutoff,]
	
	# Make an MA-plot displaying the significant reads
	bitmap(file="ma-plot-significant-counts.png", width=w/72, height=h/72)	
	significant_indices <- rownames (significant_results)
	plotSmear(dge_list, de.tags = significant_indices, main = "MA plot for significantly\ndifferentially expressed sequence tags")
	abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)
	dev.off()
}

##############################################
# Anlysis using moderated tagwise dispersion #
##############################################

if (dispersion_estimate == "tagwise") {
	# Calculate the tagwise dispersion
	number_moderating_tags <- 10
	dge_list <- estimateCommonDisp(dge_list)
	dge_list <- estimateTagwiseDisp(dge_list, prior.n = number_moderating_tags)
	
	# Statistical testing
	stat_test <- exactTest(dge_list, common.disp = FALSE)
	
	# Extract results in a nice-looking table
	number_tags <- dim (dge_list$counts) [1]
	results_table <- topTags (stat_test, n=number_tags, adjust.method="BH", sort.by="p.value")
	results_table <- results_table$table
	
	# Extract the significant tags based on adjusted p-value cutoff
	cutoff <- 0.05
	significant_results <- results_table[results_table$FDR<cutoff,]
	
	# Make an MA-plot displaying the significant reads
	bitmap(file="ma-plot-significant-counts.png", width=w/72, height=h/72)	
	significant_indices <- rownames (significant_results)
	plotSmear(dge_list, de.tags = significant_indices, main = "MA plot for significantly\ndifferentially expressed sequence tags")
	abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)
	dev.off()
}

# Create a tbale with the original counts per sample together with the statistical tests results
# ready for output in Chipster
output_table <- data.frame (dat[significant_indices,], significant_results)

# Output the table
write.table(output_table, file="de-list.tsv", sep="\t", row.names=T, col.names=T, quote=F)
