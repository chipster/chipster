# TOOL ngs-find-peaks-macs-two.R: "Find peaks using MACS, treatment vs. control" (This tool will search for statistically significantly enriched genomic regions in sequencing data from a ChIP-seq experiment. The analysis is performed on one or more treatment samples relative to one or more control samples.)
# INPUT treatment.bam: "Treatment data file" TYPE GENERIC 
# INPUT control.bam: "Control data file" TYPE GENERIC 
# OUTPUT positive-peaks.tsv: "True enriched peaks" 
# OUTPUT positive-peaks.bed: "True enriched peaks in a format compatible with the Genome Browser"
# OUTPUT OPTIONAL model-plot.pdf: "A plot of the fitted peak model" 
# OUTPUT OPTIONAL negative-peaks.tsv: "The false enriched peaks" 
# OUTPUT analysis-log.txt: "Summary of analysis settings and run" 
# PARAMETER file.format: "File format" TYPE [ELAND, BAM, BED] DEFAULT BAM (The format of the input files.)
# PARAMETER species: Genome TYPE [human, mouse, rat] DEFAULT human (the species of the samples.)
# PARAMETER version: "MACS version" TYPE [1, 2] DEFAULT 1 (Determines if analysis is done using MACS1 or MACS2)
# PARAMETER OPTIONAL read.length: "Read length" TYPE INTEGER FROM 0 TO 200 DEFAULT 0 (The length in nucleotides of the sequence reads. Option 0 envokes the default behaviour in which read length is auto-detected)
# PARAMETER OPTIONAL band.with: "Bandwidth" TYPE INTEGER FROM 1 TO 1000 DEFAULT 200 (The scanning window size, typically half the average fragment size of the DNA)
# PARAMETER OPTIONAL p.value.threshold: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.00001 (The cutoff for statistical significance. Since the p-values are not adjusted to account for multiple testing correction, the cutoff needs to be substantially more conservative than what is usually applied.)
# PARAMETER OPTIONAL build.model: "Peak model" TYPE [yes, no] DEFAULT yes (If enabled, a peak model is built from the data. Disabling model building means the shiftsize has to be guessed. In Chipster the shift size is set to half the bandwidth.)
# PARAMETER OPTIONAL keep.dup: "Keep duplicates" TYPE [auto, all, 1] DEFAULT auto (Procedure used to handle duplicate tags. If auto, MACS computes the maximum tags at the exact same location based on binomal distribution using 1e-5 p-value cutoff. All option keeps all tags, while 1 keeps only one tag per site.)
# PARAMETER OPTIONAL m.fold.upper: "Upper M-fold cutoff" TYPE INTEGER FROM 1 TO 100 DEFAULT 30 (Sets the cutoff used to determine peak regions for model building. A too high value may result in not enough peaks being identified for building the model. Notice that if the peak model is disabled this parameter has no effect.)
# PARAMETER OPTIONAL m.fold.lower: "Lower M-fold cutoff" TYPE INTEGER FROM 1 TO 100 DEFAULT 10 (Sets the cutoff used to determine peak regions for model building. A too low value may result in the inclusion of many false peaks being used for building the model. Notice that if the peak model is disabled this parameter has no effect.)

# Tool that searches for genomic regions that are significantly enriched in sequence reads from a
# ChIP-chip or ChIP-seq experiment. The tool uses the MACS algorithm for identifying the enriched
# regions and uses a p-value cutoff for determining the statistical significance. This version of 
# the tool works for experiments with two experiment condition, like treatment and control.

# 26.05.2010 MG, Created
# 01.12.2012 MG, Modified to take BAM files as input. Modified to run version 1.4 of MACS.
# 08.03.2011 MG, Modified to disable wiggle output.
# 05.04.2014 MK, Polished. Added MACS2

# MK 05.05.2014: Instead of commeting code, use version control to store old code that could be useful in future
# See version control for code that could be used if reading the experiment setup from the phenodata file, like is done for 
# microarray data. The code allows multiple samples per treatment group and will automatically merge all samples into a 
# single file per treatment group.      

# MACS settings
if(version == 1) {
	macs.binary <- file.path(chipster.tools.path, "macs", "macs142")
} else {
	macs.binary <- file.path(chipster.tools.path, "macs", "macs2")
}

# Set up approximate mappable genome size depending on species
if (species == "human") {
	genome.size <- as.character(3.5e+9*0.978*.9)
}
if (species == "mouse") {
	genome.size <- as.character(3.25e+9*0.978*.9)
}
if (species == "rat") {
	genome.size <- as.character(3.05e+9*0.978*.9)
}

if (read.length == 0) {
	read.length = FALSE
}

# Set up some parameters in case building peak model is disabled
if (build.model == "no") {
	no.model <- TRUE
	shift.size <- band.with / 2
}

# Set up some parameters in case building peak model is enabled
if (build.model == "yes") {
	no.model <- FALSE
}

# Set up the m-fold limits
mfold.limits <- paste (as.character(m.fold.lower),",",as.character(m.fold.upper), sep="")

# Define function for running MACS
runMACS <- function(..., logFile="/dev/null") {
	
	# Parameter values (character vector)
	params <- c(...) ## Nice :)
	if (!is.null(params)) {
		
		# Flags
		flags <- paste("--", names(params), sep="")
		
		# Remove the parameter names
		names(params) <- NULL
		
		# Switches that are on (value is "TRUE" or "T")
		switchOnParams <- NULL
		switchOn <- which(params == "T" | params == "TRUE")
		if (!identical(switchOn, integer(0))) {
			switchOnParams <- flags[switchOn]
			params <- params[-switchOn]
			flags <- flags[-switchOn]
		}
		
		# Switches that are off are ignored (value is "FALSE" or "F")        
		switchOff <- which(params == "F" | params == "FALSE")
		if (!identical(switchOff, integer(0))) {
			params <- params[-switchOff]
			flags <- flags[-switchOff]
		}
		
		# Command
		command <- paste(macs.binary, paste(flags, params, collapse=" "))
		if (!is.null(switchOnParams)) {
			switchOnParams <-  paste(switchOnParams, collapse=" ")
			command <- paste(command, switchOnParams)
		}
		
		return(paste(command, "2>", logFile))		
	}
}

# Run MACS with specified parameters for the data set
if (build.model == "no") {
	command <- runMACS(treatment="treatment.bam", 
						control="control.bam", 
						name="results", 
						format = file.format,
						bw=band.with,
						pvalue=p.value.threshold,
						mfold=mfold.limits,
						tsize=read.length,
						gsize=genome.size,
						"keep-dup"=keep.dup,
						verbose=3, 
						logFile="results.log", 
						nomodel=no.model,
						shiftsize=shift.size,
						help=FALSE, 
						version=FALSE)
}
if (build.model == "yes") {
	command <- runMACS(treatment="treatment.bam", 
						control="control.bam", 
						name="results", 
						format = file.format,
						bw=band.with,
						pvalue=p.value.threshold,
						mfold=mfold.limits,
						tsize=read.length,
						gsize=genome.size,
						"keep-dup"=keep.dup,
						verbose=3, 
						logFile="results.log", 
						nomodel=no.model,
						help=FALSE, 
						version=FALSE)
}


pypath <- ""
system.output <- system(gsub("^;", "", paste(pypath, command, sep=";")))

if (system.output != 0) {
	stop("CHIPSTER-NOTE: Building the peak model failed. Retry by lowering the m-fold value.") 
}
if (length(grep ("results_model.r",dir())) == 0) {
	stop("CHIPSTER-NOTE: Building the peak model failed. Retry by lowering the m.fold.lower value or rerun with model builing turned off.") 
}

# Read in and parse the results

## If final == FALSE, parse the bootsrap results. Else parse the final results
## This feature is not yet implemented.
## <name> is the experiment name and it must be same than in the runMACS
## 10xlog10(pvalue) values are sorted decreasingly"

parseMACSResults <- function(name, final=FALSE){
	if( final ) {
		output <- read.table(file=name, skip=0, header=TRUE, stringsAsFactors=FALSE)
		colnames(output)[7] <- "neg10xlog10pvalue"
		## Sort the result according to the -10xlog10(pvalue)
		output <- output[ order(output[,5], decreasing=TRUE), ]
		return(output)
	}
}

## Read in the results for the TRUE peaks
results_TRUE <- parseMACSResults (name="results_peaks.xls",final=TRUE)
results_TRUE$chr <- sub(".fa", "", results_TRUE$chr)
results_TRUE$chr <- sub("chr", "", results_TRUE$chr)
results_TRUE_ordered <- results_TRUE[order(results_TRUE$chr, results_TRUE$start),]
write.table(results_TRUE_ordered, file="positive-peaks.tsv", sep="\t", quote=FALSE, row.names=FALSE)

if (version == 1) {
	# MACS1 computes empirical FDR by swapping ChIP and control. These are stored in "negative_peaks" file 
	results_FALSE <- parseMACSResults (name="results_negative_peaks.xls", final=TRUE)
	results_FALSE$chr <- sub(".fa", "", results_FALSE$chr)
	results_FALSE$chr <- sub("chr", "", results_FALSE$chr)
	results_FALSE_ordered <- results_FALSE[order(results_FALSE$chr, results_FALSE$start),]
	write.table(results_FALSE_ordered, file="negative-peaks.tsv", sep="\t", quote=FALSE, row.names=FALSE)
}

# Convert the name of some files to make it compatible with chipster output
system("mv results.log analysis-log.txt")

# Sort the BED
source(file.path(chipster.common.path, "bed-utils.R"))
if (file.exists("results_peaks.bed")){
	bed <- read.table(file="results_peaks.bed", skip=0, sep="\t")
	colnames(bed)[1:2] <- c("chr", "start")
	bed <- sort.bed(bed)
	write.table(bed, file="positive-peaks.bed", sep="\t", row.names=F, col.names=F, quote=F)
}

# Source the R code for plotting the MACS model and convert the PDF file to PNG
if (build.model == "yes") {
	source("results_model.r")
	system("mv results_model.pdf model-plot.pdf")
}
