# TOOL "Statistics" / ngs-find-motifs-jaspar.R: "Find common motifs and match to Jaspar" (This tool scans a set of genomic regions for consensus sequence motifs,
# calculates the alignment score against transcription factors in the Jaspar database and finds the 10 highest ranking for each motif.)
# INPUT results.tsv: "Results data file" TYPE GENERIC
# OUTPUT motif-analysis-summary.txt: "A lot of analysis information collected in one single file"
# OUTPUT logo-plot-(...).png: "Logo plots for each consensus motif"
# OUTPUT test_plot_(...).png: "testing optional plots"
# OUTPUT testplot.png: "testing optional plots"
# PARAMETER p.value.cutoff: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.0002 (This parameter controls the false positive rate when searching for consensus sequence motifs. Lower the value for increased stringency.)
# PARAMETER e.value.cutoff: "E-value cutoff" TYPE DECIMAL FROM 0 TO 100 DEFAULT 0.01 (This parameter controls the alignment stringency, where a lower value means better alignment.)
# PARAMETER genome: "Genome" TYPE [BSgenome.Hsapiens.UCSC.hg17, BSgenome.Hsapiens.UCSC.hg18, BSgenome.Hsapiens.UCSC.hg19, BSgenome.Mmusculus.UCSC.mm8, BSgenome.Mmusculus.UCSC.mm9, BSgenome.Rnorvegicus.UCSC.rn4] DEFAULT BSgenome.Hsapiens.UCSC.hg18 (The genome and version used when aligning the sequences.)

#####################################################
#                                                   #
# MG, 26.5.2010                                     #
#                                                   #
# Development version                               #
#                                                   #
# Tool that searches for consensus motifs in a set  #
# of genomic regions, typically a result from a     #
# peak finding analysis for ChIP-seq data,          #
# calculaates alignemt scores to TF from the JASPAR #
# database and collects results for the 10 highest  #
# TF for each motif.                                #
#                                                   #
#####################################################


# Load the required libraries
library(MotIV)
library(rGADEM)
library(package=genome, character.only=TRUE)

# Covert genome name to fit requirements from rGADEM package
#if (genome == "BSgenome.Hsapiens.UCSC.hg17" | genome == "BSgenome.Hsapiens.UCSC.hg18" | genome == "BSgenome.Hsapiens.UCSC.hg19") {
#		genome <- "Hsapiens"
#}
#if (genome == "BSgenome.Mmusculus.UCSC.mm8" | genome == "BSgenome.Mmusculus.UCSC.mm9") {
#	genome <- "Mmusculus"
#}
#if (genome == "BSgenome.Rnorvegicus.UCSC.rn4") {
#	genome <- "Rnorvegicus"
#}

# Read in data and convert to BED format
results_file <- read.table (file="results.tsv", sep="\t", header=T)
results_bed <- results_file[,1:3]
results_bed[,1] <- paste("chr",results_bed[,1], sep="")

# Convert to Ranged data
results_ranged <- IRanges(start=results_bed[,2], end=results_bed[,3])
results_sequences <- RangedData(results_ranged, space=results_bed[,1])

# Perform unseeded GADEM analysis with default values for specific genome
if (genome == "BSgenome.Hsapiens.UCSC.hg17" | genome == "BSgenome.Hsapiens.UCSC.hg18" | genome == "BSgenome.Hsapiens.UCSC.hg19") {
	results_gadem <- GADEM(
			results_sequences,
			verbose=1,
			genome=Hsapiens,
			pValue=p.value.cutoff,
			eValue=e.value.cutoff)
}
if (genome == "BSgenome.Mmusculus.UCSC.mm8" | genome == "BSgenome.Mmusculus.UCSC.mm9") {
	results_gadem <- GADEM(
			results_sequences,
			verbose=1,
			genome=Mmusculus,
			pValue=p.value.cutoff,
			eValue=e.value.cutoff)
}
if (genome == "BSgenome.Rnorvegicus.UCSC.rn4") {
	results_gadem <- GADEM(
			results_sequences,
			verbose=1,
			genome=Rnorvegicus,
			pValue=p.value.cutoff,
			eValue=e.value.cutoff)
}

# Read in Jaspar database
# path_jaspar <- system.file(package="rGADEM")
# file_name <- paste(path_jaspar, "/extdata/jaspar2009.txt", sep="/")
# jaspar <- readPWMfile(file_name)
#
# Actually, as of veersion 1.1.6 of MotIV the Jaspar database and scores
# are loaded automatically whem MotIV is loaded
# the data is called jaspar and the scores jaspar.scores

# Get the motifs from the GADEM objects
results_motifs <- viewPWM(results_gadem)

# Find out how many consensus motifs were discovered
number_motifs <- length(results_motifs)

# Find our how many occurrences of each motif
number_occurrences <- nOccurrences (results_gadem)

# View the consensus sequence for each motif
consensus_sequences <- consensus (results_gadem)

# Give appropriate names to the motifs
names <- character (number_motifs)
for (count in 1:number_motifs) {
	names [count] <- paste("motif", count, sep=" ")
}
names (results_motifs) <- names

# Perform a MotIV analysis of alignment with TF:s from JASPAR, only the top 10 matches for each motif are collected
results_alignment <- motifMatch(
		inputPWM=results_motifs,
		align="SWU",
		cc="PCC",
		database=jaspar,
		DBscores=jaspar.scores,
		top=10)

# Get a list of the TF that matched any of the motifs
matching_tf_list <- viewMotifs (results_alignment)

# Get a summary of the results and settings
results_summary <- summary (results_alignment)

# Write out a results summary file
# print out a summary of the results
sink(file="motif-analysis-summary.txt")
	print("Summary of MotIV analysis", quote=FALSE)
	print("", quote=FALSE)
	summary(results_alignment)
	print("", quote=FALSE)
	print("", quote=FALSE)
	print ("Consensus sequence", quote=FALSE)
	print("", quote=FALSE)
	consensus_sequences
	print("", quote=FALSE)
	print("", quote=FALSE)
	print ("Number of occurrences of consensus sequence for each motif", quote=FALSE)
	print("", quote=FALSE)
	number_occurrences
	print("", quote=FALSE)
	print("", quote=FALSE)
	print ("Weighted matrices for each motif", quote=FALSE)
	print("", quote=FALSE)
	results_motifs
sink()

# testing optional plots
x <- rnorm (100,50,10)
y <- rnorm (100,25,5)
bitmap(file="test_plot_1.png", width=1000/72, height=1000/72)
plot(x,y)
dev.off()
bitmap(file="test_plot_2.png", width=1000/72, height=1000/72)
plot(y,x)
dev.off()
bitmap(file="testplot.png", width=1000/72, height=1000/72)
plot(y,x)
dev.off()

# Plot the logo of the motifs and corresponding TF:s
for (count in 1:number_motifs) {
	file_name <- paste("logo-plot-",count,".png", sep="")
	bitmap(file=file_name, width=1000, height=1000)
#	png(width=1000, height=1000, file=file_name)
	plot_name <- consensus(results_gadem) [count]
	plot (results_alignment[count], top=10, main=paste("Top ten TF matches for motif\n", plot_name, sep=""))
	dev.off()
}

