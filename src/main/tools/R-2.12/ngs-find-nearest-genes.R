# TOOL "Statistics" / ngs-find-nearest-genes.R: "Find the nearest genes for regions" (This tool takes set of genomic regions, such as ChIP-seq peaks, and fetches the nearest gene, exon or miRNA for each.)
# INPUT regions-list.tsv: "Table with genomic regions" TYPE GENERIC
# OUTPUT nearest-genes.tsv: "Table listing the nearest gene feature for each input region."
# PARAMETER species: "Species" TYPE [Human, Mouse, Rat] DEFAULT Human (The species of the genome to use for fetching annotationsan.)

#####################################################
#                                                   #
# MG, 27.10.2010                                    #
#                                                   #
# Development version                               #
#                                                   #
# Tool that fetches the nearest gene, exon or miRNA #
# for a set of genomic regions, such as the output  #
# from a peak detection algorithm for ChIP-seq data #
#                                                   #
#####################################################

# Possible additional parameters and output
# OUTPUT logo-plot-{...}.png: "Logo plots for each consensus motif"
# PARAMETER feature: "Genomic feature" TYPE [TSS, exon, miRNA] DEFAULT TSS (The genomic feature to look for in each region. Currently, transcription start site, exon or miRNA are supported feature types)
# PARAMETER e.value.cutoff: "E-value cutoff" TYPE DECIMAL FROM 0 TO 100 DEFAULT 0.01 (This parameter controls the alignment stringency, where a lower value means better alignment.)

# Set up testing parameters
species <- "Human"

# Load the libraries
library(ChIPpeakAnno)

# Load the annotation data
if (species == "Human") {
	data(TSS.human.NCBI36)
	annotations <- TSS.human.NCBI36
}
if (species == "Nouse") {
	data(TSS.mouse.NCBIM37)
	annotations <- TSS.mouse.NCBIM37
}
if (species == "Rat") {
	data(TSS.rat.RGSC3.4)
	annotations <- TSS.rat.RGSC3.4
}

# Read in data and convert to BED format
results_file <- read.table (file="regions-list.tsv", sep="\t", header=T)
results_bed <- results_file[,1:3]
results_bed[,1] <- paste("chr",results_bed[,1], sep="")

# Convert to Ranged data
results_ranged <- BED2RangedData(results_bed)

# Annotate the data
results_annotated = annotatePeakInBatch(results_ranged, AnnotationData = TSS.human.NCBI36)

# Convert into a table and extract the useful columns
results_table <- as.data.frame(results_annotated)
rownames(results_table) <- sub (" ", "_", results_table$name)
results_table <- results_table[,c(-5,-13,-14)]
table_header <- c("chromosome","region_start","region_end","region_width",
	"region_id","strand","ensembl_id","gene_start","gene_end","location",
	"distance")
names(results_table) <- table_header

# Write output
write.table(results_table, file="nearest-genes.tsv", sep="\t", col.names=T, row.names=T, quote=F)

# EOF