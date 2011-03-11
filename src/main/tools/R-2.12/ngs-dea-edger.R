# ANALYSIS Statistics/"Two groups tests" (Tests for comparing the mean gene expression of two groups. 
# LPE only works, if the whole data is used, i.e., the data should not be pre-filtered, if LPE is used. 
# Other than empiricalBayes might be slow, if run on unfiltered data.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT two-sample.tsv
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER test [empiricalBayes, fast-t-test, t-test, F-test, Mann-Whitney, LPE] DEFAULT empiricalBayes (Test type)
# PARAMETER p.value.adjustment.method [none, Bonferroni, Holm, Hochberg, BH, BY] DEFAULT BH (Multiple testing correction method)
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)

# TOOL "Statistics" / ngs-find-peaks-macs-one.R: "Find peaks using MACS, treatment only" (This tool will search for statistically significantly enriched genomic regions in sequencing data from a ChIP-seq experiment. The analysis is performed on one or more treatment samples alone, without taking into account control control samples.)
# INPUT treatment.bam: "Treatment data file" TYPE GENERIC
# OUTPUT positive-peaks.tsv: "True enriched peaks"
# OUTPUT OPTIONAL model-plot.png: "A plot of the fitted peak model"
# OUTPUT OPTIONAL negative-peaks.tsv: "The false enriched peaks"
# OUTPUT analysis-log.txt: "Summary of analysis settings and run"
# PARAMETER file.format: "File format" TYPE [ELAND, BAM, BED] DEFAULT BAM (The format of the input files.)
# PARAMETER species: "Genome" TYPE [human, mouse, rat] DEFAULT human (the species of the samples.)
# PARAMETER read.length: "Read length" TYPE INTEGER FROM 1 TO 200 DEFAULT 25 (The length in nucleotides of the sequence reads)
# PARAMETER band.with: "Band width" TYPE INTEGER FROM 1 TO 1000 DEFAULT 200 (The scanning window size, typically half the average fragment size of the DNA)
# PARAMETER p.value.threshold: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.00001 (The cutoff for statistical significance. Since the p-values are not adjusted to account for multiple testing correction the cutoff needs to be substantially more conservative than what is usually applied.)
# PARAMETER build.model: "Peak model" TYPE [yes, no] DEFAULT yes (If enabled, a peak model is built from the data. Disabling model building means the shiftsize has to be guessed. In Chipster the shift size is set to half the band with.)
# PARAMETER m.fold.upper: "Upper M-fold cutoff" TYPE INTEGER FROM 1 TO 100 DEFAULT 30 (Sets the cutoff used to determine peak regions for model building. A too high value may result in not enough peaks being identified for building the model. Notice that if the peak model is disabled this parameter has no effect.)
# PARAMETER m.fold.lower: "Lower M-fold cutoff" TYPE INTEGER FROM 1 TO 100 DEFAULT 10 (Sets the cutoff used to determine peak regions for model building. A too low value may result in the inclusion of many false peaks being used for building the model. Notice that if the peak model is disabled this parameter has no effect.)


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
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Test needs a parameter "groups" that specifies the grouping of the samples
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
groups<-phenodata[,pmatch(column,colnames(phenodata))]

# Sanity checks
if(length(unique(groups))==1 | length(unique(groups))>=3) {
   stop("You need to have exactly two groups to run this analysis")
}



# Read in the data from one sample
setwd("C:/Chipster testing/edgeR")
data <- read.table(file="peaksHigherThan10", sep="\t", header=F, quote="")

# Remove rows pertaining to "*" chromosome
sample_1 <- data[grep ("*",as.character (data[,2]), useBytes=TRUE, fixed=TRUE, ivert=TRUE),]

# Create 5 artificial samples by randomly sampling values from sample 1
sample_2 <- sample_1[sample(1:dim(sample_1)[1], size=dim(sample_1)[1], replace=FALSE),4]
sample_3 <- sample_1[sample(1:dim(sample_1)[1], size=dim(sample_1)[1], replace=FALSE),4]
sample_4 <- sample_1[sample(1:dim(sample_1)[1], size=dim(sample_1)[1], replace=FALSE),4]
sample_5 <- sample_1[sample(1:dim(sample_1)[1], size=dim(sample_1)[1], replace=FALSE),4]
sample_6 <- sample_1[sample(1:dim(sample_1)[1], size=dim(sample_1)[1], replace=FALSE),4]

# Bind the samples together in the same datatable
data_2 <- cbind(sample_1, sample_2)
data_2 <- cbind(data_2, sample_3)
data_2 <- cbind(data_2, sample_4)
data_2 <- cbind(data_2, sample_5)
data_2 <- cbind(data_2, sample_6)

# Fix column names
colnames(data_2)[c(1,2,3,4)] <- c("seq_id","chr","start","sample_1")

# Create a "phenodata" type file to describe the experiment
pheno_data <- c(rep("experiment",3),rep("control",3))

# Create a DGEList
# Notice that Library size is calculated from column totals
dge_list <- DGEList (count=data_2[,4:9], group=pheno_data)

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
write.table(output_table, digits=6)), file="two-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)

