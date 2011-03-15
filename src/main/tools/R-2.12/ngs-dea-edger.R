# TOOL "Statistics" / ngs-dea-edger.R: "Differential expression analysis using edgeR" (This tool will perform an analysis for differentially expressed sequences using the R implementation of the edge algorithm.)
# INPUT data.tsv: "Input data table" TYPE GENE_EXPRS
# INPUT phenodata.tsv "Experiment parameters" TYPE GENERIC
# OUTPUT de-list.tsv "List of differentially expressed sequences" TYPE GENE_EXPR
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER normalization_method: "Normalization method" TYPE [none, library_size, sequence_specific] DEFAULT sequence_specific (If enabled, a normalization factor based on the total number of sequence counts or sequence distribution is performed to reduce the effect from sequencing biases.)
# PARAMETER dispersion_estimate: "Dispersion estimate" TYPE [common, moderated] DEFAULT moderated (The dispersion of counts for any given sequence can either be estimated based on the actual counts in the sample data set or be moderated across a selection of sequences with similar count numbers. The latter option, which is set by default, typically yields higher sensitivity and specificity.)
# PARAMETER p_value_adjustment_method "Multiple testing correction" TYPE [none, Bonferroni, Holm, Hochberg, BH, BY] DEFAULT BH (Multiple testing correction method.)
# PARAMETER p_value_threshold: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.00001 (The cutoff for statistical significance.)


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

# Loads the normalized data
file <- c("normalized.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls <- dat[,grep("flag", names(dat))]
dat2 <- dat[,grep("chip", names(dat))]

# Test needs a parameter "groups" that specifies the grouping of the samples
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
groups <- phenodata[,pmatch(column,colnames(phenodata))]

# Sanity checks
if(length(unique(groups))==1 | length(unique(groups))>=3) {
   stop("You need to have exactly two groups to run this analysis")
}

# Create a "phenodata" type file to describe the experiment
pheno_data <- c(rep("experiment",3),rep("control",3))

# Create a DGEList
# Notice that Library size is calculated from column totals
dge_list <- DGEList (count=data_2[,4:9], group=groups)

# Check lib size totals
dge_list$samples

# Calculate normalization factors
dge_list <- calcNormFactors(dge_list) 

# MA-plot comparison before and after normalization
par(mfrow = c(1, 2))
maPlot(dge_list$counts[, 1], dge_list$counts[, 2], normalize = TRUE, pch = 19,
		cex = 0.4, ylim = c(-8, 8))
grid(col = "blue")
title("Before normalization")
abline(h = log2(dge_list$samples$norm.factors[2]/dge_list$samples$norm.factors[1]),
		col = "red", lwd = 4)
eff.libsize <- dge_list$samples$lib.size * dge_list$samples$norm.factors
maPlot(dge_list$counts[, 1]/eff.libsize[1], dge_list$counts[, 2]/eff.libsize[2],
		normalize = FALSE, pch = 19, cex = 0.4, ylim = c(-8, 8))
grid(col = "blue")
title("After normalization")

# Produce MDS plot of normazied data
sample_colors <-  ifelse (dge_list$samples$group=="experiment", 1, 2)
plotMDS.dge(dge_list, main="MDS Plot", xlim=c(-2,1), col=sample_colors)
legend(x="topleft", legend = c("experiment","control"),col=c(1,2), pch=19)


###################################
# Anlysis using common dispersion #
###################################

# Calculate common dispersion
dge_list <- estimateCommonDisp(dge_list)

# Statistical testing
stat_test <- exactTest(dge_list)

# Extract results in a nice-looking table
number_tags <- dim (dge_list$counts) [1]
results_table <- topTags (stat_test, n=number_tags, adjust.method="BH", sort.by="p.value")
results_table <- results_table$table

# Extract the significant tags based on adjusted p-value cutoff
cutoff <- 0.05
significant_results <- results_table[results_table$FDR<cutoff,]

# Make an MA-plot displaying the significant reads
significant_indices <- rownames (significant_results)
plotSmear(dge_list, de.tags = significant_indices, main = "MA plot using common dispersion")
abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)


##############################################
# Anlysis using moderated tagwise dispersion #
##############################################

# Calculate the tagwise dispersion
number_moderating_tags <- 10
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
significant_indices <- rownames (significant_results)
plotSmear(dge_list, de.tags = significant_indices, main = "MA plot using tagwise dispersion")
abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)

# Create a tbale with the original counts per sample together with the statistical tests results
# ready for output in Chipster
output_table <- data.frame (data_2[significant_indices,], significant_results)

# Output the table
write.table(output_table, file="two-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)

