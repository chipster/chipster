# TOOL ngs-find-peaks-macs-two.R: "Find peaks using MACS, treatment vs. control" (Detects statistically significantly enriched genomic regions in ChIP-seq data using a control sample. If you have several samples, you need to merge them first to one ChIP file and one control file. BAM files can be merged with the Utilities tool \"Merge BAM\".)
# INPUT treatment.bam: "Treatment data file" TYPE GENERIC 
# INPUT control.bam: "Control data file" TYPE GENERIC 
# OUTPUT positive-peaks.tsv: "True enriched peaks" 
# OUTPUT positive-peaks.bed: "True enriched peaks in a format compatible with the Genome Browser"
# OUTPUT OPTIONAL model-plot.pdf: "A plot of the fitted peak model" 
# OUTPUT OPTIONAL negative-peaks.tsv: "The false enriched peaks" 
# OUTPUT analysis-log.txt: "Summary of analysis settings and run" 
# PARAMETER file.format: "Input file format" TYPE [ELAND, BAM, BED] DEFAULT BAM (The format of the input files.)
# PARAMETER precalculated.size: "Mappable genome size" TYPE [2.7e9: "human hg18 (2.7e9\)", 2.72e9: "human hg19 (2.72e9\)", 1.87e9: "mouse mm9 (1.87e9\)", 1.89e9: "mouse mm10 (1.89e9\)", 2.32e9: "rat rn5 (2.32e9\)", user_specified: "User specified"] DEFAULT 2.72e9 (Mappable genome size. You can use one of the precalculated ones or choose User specified and provide the size in the field below.)
# PARAMETER OPTIONAL userspecifed.size: "User specified mappable genome size" TYPE STRING (You can also use scientific notation, e.g. 1.23e9 . Remember to select User specified as Mappable genome size.)
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
# 10.07.2014 AMS, Updated genome sizes, added parameter userspecified.size
# 10.09.2014 EK, Removed MACS2 to a separate tool, fixed the bug in disabled model building, polished the script and output

# MK 05.05.2014: 
# See version control for code that could be used if reading the experiment setup from the phenodata file, like is done for microarray data.      


# MACS binary
macs.binary <- file.path(chipster.tools.path, "macs", "macs14")

# Use user-specified genome size if given
if (precalculated.size == "user_specified") {
	if (nchar(userspecifed.size) < 1){
		stop(paste('CHIPSTER-NOTE: ', "You need to provide a value for mappable genome size or select one of the precalculated values."))
	}
	genome.size <- userspecifed.size
}else{
	genome.size <- precalculated.size
}

# If read length is left to zero, it should be estimated from data
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
						verbose=2, 
						logFile="results.log", 
						nomodel=no.model,
						shiftsize=shift.size)
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
						verbose=2, 
						logFile="results.log", 
						nomodel=no.model)
}


pypath <- ""
system.output <- system(gsub("^;", "", paste(pypath, command, sep=";")))

if (system.output != 0) {
	stop("CHIPSTER-NOTE: Building the peak model failed. Retry by lowering the m-fold value.") 
}
#if (length(grep ("results_model.r",dir())) == 0) {
#	stop("CHIPSTER-NOTE: Building the peak model failed. Retry by lowering the m.fold.lower value or rerun with model building turned off.") 
#}

# Read in and parse the results, sort 
output <- read.table(file="results_peaks.xls", skip=0, header=TRUE, stringsAsFactors=FALSE)
colnames(output)[7] <- "neg10xlog10pvalue"
colnames(output)[9] <- "FDR_percentage"
# output <- output[ order(output[,7], decreasing=TRUE), ]
output <- output[order(output$chr, output$start),]
write.table(output, file="positive-peaks.tsv", sep="\t", quote=FALSE, row.names=FALSE)

negoutput <- read.table(file="results_negative_peaks.xls", skip=0, header=TRUE, stringsAsFactors=FALSE)
colnames(negoutput)[7] <- "neg10xlog10pvalue"
# negoutput <- negoutput[ order(negoutput[,7], decreasing=TRUE), ]
output <- output[order(output$chr, output$start),]
write.table(negoutput, file="negative-peaks.tsv", sep="\t", quote=FALSE, row.names=FALSE)

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
