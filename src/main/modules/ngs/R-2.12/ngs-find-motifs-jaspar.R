# TOOL ngs-find-motifs-jaspar.R: "Find common motifs and match to JASPAR" (This tool scans a set of genomic regions for consensus sequence motifs,
# calculates the alignment score against transcription factors in the JASPAR database, and reports the 10 highest ranking for each motif.)
# INPUT results.tsv: "Results data file" TYPE GENERIC 
# OUTPUT motif-analysis-summary.txt: "A lot of analysis information collected in one single file" 
# OUTPUT logo-plot-{...}.pdf: "Logo plots for each consensus motif" 
# PARAMETER p.value.cutoff: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.0002 (This parameter controls the false positive rate when searching for consensus sequence motifs. Lower the value for increased stringency.)
# PARAMETER e.value.cutoff: "E-value cutoff" TYPE DECIMAL FROM 0 TO 100 DEFAULT 0.01 (This parameter controls the alignment stringency, where a lower value means better alignment.)
# PARAMETER genome: Genome TYPE [BSgenome.Hsapiens.UCSC.hg17: hg17, BSgenome.Hsapiens.UCSC.hg18: hg18, BSgenome.Hsapiens.UCSC.hg19: hg19, BSgenome.Mmusculus.UCSC.mm8: mm8, BSgenome.Mmusculus.UCSC.mm9: mm9, BSgenome.Rnorvegicus.UCSC.rn4: rn4] DEFAULT BSgenome.Hsapiens.UCSC.hg18 (The genome and version used when aligning the sequences.)
# PARAMETER number.best.matches: "Number of matches" TYPE INTEGER FROM 1 TO 20 DEFAULT 10 (The number of best matching transcription factors for each consensus sequence found. This affects both the textual summary output and the LOGO plots.)

# MG, 26.5.2010
# MG, 6.10.2011, added parameter to control number of best matches per TF amd updated to changes in R-2.12.1 and a

# Load the required libraries
library(MotIV)
library(rGADEM)
library(package=genome, character.only=TRUE)

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
# but first make sure that there were some motifs found
# motifs_found <- try (results_motifs <- viewPWM(results_gadem))
#if (class(motifs_found) == "try-error") {
#	stop("CHIPSTER-NOTE: No common motifs were found among the query sequences! Retry with less stringent parameter settings or provide a longer list of query sequences.")
#}
if (length(nOccurrences(results_gadem)) < 1) {
	stop("CHIPSTER-NOTE: No common motifs were found among the query sequences! Retry with less stringent parameter settings or provide a longer list of query sequences.")
}

# Extract the PWM from the motif analysis
pwm_list <- getPWM(results_gadem)

# Find out how many consensus motifs were discovered
number_motifs <- length (nOccurrences (results_gadem))

# Find our how many occurrences of each motif
number_occurrences <- nOccurrences (results_gadem)

# View the consensus sequence for each motif
consensus_sequences <- consensus (results_gadem)

# Give appropriate names to the motifs
names <- character (number_motifs)
for (count in 1:number_motifs) {
	names [count] <- paste("motif", count, sep=" ")
}
names (pwm_list) <- names

# Perform a MotIV analysis of alignment with TF:s from JASPAR, only the top 10 matches for each motif are collected
results_alignment <- motifMatch(
		inputPWM=pwm_list,
		align="SWU",
		cc="PCC",
		database=jaspar,
		DBscores=jaspar.scores,
		top=number.best.matches)

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
print ("Alignments of motifs to known TFBS", quote=FALSE)
print("", quote=FALSE)
viewAlignments(results_alignment)
print("", quote=FALSE)
print("", quote=FALSE)
print ("Weighted matrices for each motif", quote=FALSE)
print("", quote=FALSE)
pwm_list
sink()

# Plot the logo of the motifs and corresponding TF:s
# Adjust y-scale to number of best matches to be plotted
if (number.best.matches < 11) { 
	for (count in 1:number_motifs) {
		file_name <- paste("logo-plot-",count,".pdf", sep="")
		pdf(file=file_name, height=1400/72, width=600/72)
		plot_name <- consensus(results_gadem) [count]
		plot (results_alignment[count], top=number.best.matches, ncol=1, main=paste("Top TF matches for motif\n", plot_name, sep=""))
		dev.off()
	}
} else {
	for (count in 1:number_motifs) {
		file_name <- paste("logo-plot-",count,".pdf", sep="")
		pdf(file=file_name, height=2000/72, width=600/72)
		plot_name <- consensus(results_gadem) [count]
		plot (results_alignment[count], top=number.best.matches, ncol=1, main=paste("Top TF matches for motif\n", plot_name, sep=""))
		dev.off()
	}
}

# EOF
