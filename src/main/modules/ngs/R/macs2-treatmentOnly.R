# TOOL macs2-treatmentOnly.R: "Find peaks using MACS2, treatment only" (Detects statistically significantly enriched genomic regions in ChIP-seq data. The analysis does not use a control sample. If you have several ChIP samples,you need to merge them first to one file. BAM files can be merged with the Utilities tool \"Merge BAM\".)
# INPUT treatment.bam: "Treatment data file" TYPE GENERIC 
# OUTPUT macs2-peaks.tsv
# OUTPUT macs2-log.txt
# OUTPUT macs2-peaks.bed
# OUTPUT macs2-summits.bed
# OUTPUT OPTIONAL macs2_narrowpeak.bed
# OUTPUT OPTIONAL macs2_broad_peaks.bed
# OUTPUT OPTIONAL macs2_model.pdf  
# PARAMETER file.format: "Input file format" TYPE [ELAND, BAM, BED] DEFAULT BAM (The format of the input files.)
# PARAMETER precalculated.size: "Mappable genome size" TYPE [2.7e9: "human hg18 (2.7e9\)", 2.72e9: "human hg19 (2.72e9\)", 1.87e9: "mouse mm9 (1.87e9\)", 1.89e9: "mouse mm10 (1.89e9\)", 2.32e9: "rat rn5 (2.32e9\)", user_specified: "User specified"] DEFAULT 2.72e9 (Mappable genome size. You can use one of the precalculated ones or choose User specified and provide the size in the field below.)
# PARAMETER OPTIONAL userspecified.size: "User specified mappable genome size" TYPE STRING (You can also use scientific notation, e.g. 1.23e9 . Remember to select User specified as Mappable genome size.)
# PARAMETER OPTIONAL q.value.threshold: "q-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.01 (The minimum FDR for peak detection.)
# PARAMETER OPTIONAL read.length: "Read length" TYPE INTEGER FROM 0 TO 200 DEFAULT 0 (The length in nucleotides of the sequence reads. Option 0 envokes the default behaviour in which read length is auto-detected)
# PARAMETER OPTIONAL keep.dup: "Keep duplicates" TYPE [auto, all, 1] DEFAULT auto (Procedure used to handle duplicate tags. If auto, MACS computes the maximum tags at the exact same location based on binomal distribution using 1e-5 p-value cutoff. All option keeps all tags, while 1 keeps only one tag per site.)
# PARAMETER OPTIONAL build.model: "Build peak model" TYPE [yes, no] DEFAULT yes (If enabled, a peak model is built from the data. Disabling model building means the shiftsize has to be guessed and set with the parameter.)
# PARAMETER OPTIONAL bandwidth: "Bandwidth" TYPE INTEGER FROM 1 TO 1000 DEFAULT 300 (The scanning window size, typically half the average fragment size of the DNA.)
# PARAMETER OPTIONAL m.fold.upper: "Upper M-fold cutoff" TYPE INTEGER FROM 1 TO 100 DEFAULT 30 (Sets the cutoff used to determine peak regions for model building. A too high value may result in not enough peaks being identified for building the model. Notice that if the peak model is disabled this parameter has no effect.)
# PARAMETER OPTIONAL m.fold.lower: "Lower M-fold cutoff" TYPE INTEGER FROM 1 TO 100 DEFAULT 10 (Sets the cutoff used to determine peak regions for model building. A too low value may result in the inclusion of many false peaks being used for building the model. Notice that if the peak model is disabled this parameter has no effect.)                     
# PARAMETER OPTIONAL shift.size: "Shift size" TYPE INTEGER FROM 1 TO 1000 DEFAULT 100 (When model building has been switched off or when it fails, MACS will use this value as half of the fragment size to shift and extend reads.)
# PARAMETER OPTIONAL broad: "Call broad peaks" TYPE [yes, no] DEFAULT no (Call broad peaks by linking nearby highly enriched regions.)

# 26.05.2010 MG, Created
# 24.11.2010 MG, Modified to take BAM files as input. Modified to run version 1.4 of MACS.
# 08.03.2011 MG, Modified to disable wiggle output.
# 05.05.2014 MK, Script polishing. Added MACS2
# 10.07.2014 AMS, Updated genome sizes, added parameter userspecified.size
# 12.09.2014 EK, Made a separate script for MACS2 in order to cope with new parameters, added the broad option and outputs, polished the script and output

# MACS settings
macs.binary <- file.path(chipster.tools.path, "macs", "macs2")


# Use user-specified genome size if given
if (precalculated.size == "user_specified") {
	if (nchar(userspecified.size) < 1){
		stop(paste('CHIPSTER-NOTE: ', "You need to provide a value for mappable genome size or select one of the precalculated values."))
	}
	genome.size <- userspecified.size
}else{
	genome.size <- precalculated.size
}

# If read length is left to zero, it should be estimated from data
if (read.length == 0) {
	read.length = FALSE
}

# Set up some parameters
if (build.model == "no") {
	no.model <- TRUE
}
if (build.model == "yes") {
	no.model <- FALSE
}
if (broad == "yes") {
	call.broad <- TRUE
}
if (broad == "no") {
	call.broad <- FALSE
}

# Set up the m-fold limits
mfold.limits <- paste (as.character(m.fold.lower),",",as.character(m.fold.upper), sep="")

# MK 05.05.2014: Instead of commenting code, use version control to store old code that could be useful in future
# See version control for code that could be used if reading the experiment setup from the phenodata file, like is done for 
# microarray data. The code allows multiple samples per treatment group and will automatically merge all samples into a 
# single file per treatment group.      

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

		# Run macs. Macs writes its output to stderr (stream number 2)
		# &> redirects both stderr and stdout
		# Iterates through mfold values to find low enough that works
		# if (build.model == "yes" & adjust.mfold == "yes") {
		#	for (mfold in list(32, 24, 16, 8)) {
		#		system.output <- system(paste(command, paste("--mfold=", mfold, sep=""), "2>", logFile))
		#		if (system.output == 0) {
		#			break; # was succesfull, don't lower mfold value any more
		#		}
		#	}
		#}

		return(paste(command, "2>", logFile))
	}
}

# Run MACS with specified parameters for the data set
if (build.model == "no") {
	command <- runMACS(treatment="treatment.bam", 
						name="macs2", 
						format = file.format,
						bw=bandwidth,
						qvalue=q.value.threshold,
						mfold=mfold.limits,
						tsize=read.length,
						gsize=genome.size,
						verbose=2, 
						"keep-dup"=keep.dup,
						logFile="macs2-log.txt", 
						nomodel=no.model,
						shiftsize=shift.size,
						broad=call.broad)
}
if (build.model == "yes") {
	command <- runMACS(treatment="treatment.bam", 
						name="macs2", 
						format = file.format,
						bw=bandwidth,
						qvalue=q.value.threshold,
						mfold=mfold.limits,
						tsize=read.length,
						gsize=genome.size,
						"keep-dup"=keep.dup,
						verbose=2, 
						logFile="macs2-log.txt", 
						nomodel=no.model,
						"auto-bimodal"=TRUE,
						shiftsize=shift.size,
						broad=call.broad)
}

pypath <- ""
system.output <- system(gsub("^;", "", paste(pypath, command, sep=";")))

print(command)
if (system.output != 0) {
	stop("CHIPSTER-NOTE: Building the peak model failed. Retry by lowering the m-fold value.") 
}

# Read in and parse the results
output <- read.table(file="macs2_peaks.xls", skip=0, header=TRUE, stringsAsFactors=FALSE)
colnames(output)[7] <- "neglog10pvalue"
colnames(output)[9] <- "neglog10qvalue"
output <- output[order(output$chr, output$start),]
write.table(output, file="macs2-peaks.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Sort the peaks BED
source(file.path(chipster.common.path, "bed-utils.R"))
if (file.exists("macs2_peaks.bed")){
	bed <- read.table(file="macs2_peaks.bed", skip=0, sep="\t")
	colnames(bed)[1:2] <- c("chr", "start")
	bed <- sort.bed(bed)
	write.table(bed, file="macs2-peaks.bed", sep="\t", row.names=F, col.names=F, quote=F)
}

# Sort the summit BED
source(file.path(chipster.common.path, "bed-utils.R"))
if (file.exists("macs2_summits.bed")){
	bed <- read.table(file="macs2_summits.bed", skip=0, sep="\t")
	colnames(bed)[1:2] <- c("chr", "start")
	bed <- sort.bed(bed)
	write.table(bed, file="macs2-summits.bed", sep="\t", row.names=F, col.names=F, quote=F)
}

# Add BED extension to the narrow peak format file
system ("mv macs2_peaks.encodePeak macs2_narrowpeak.bed")


# Source the R code for plotting the MACS model and convert the PDF file to PNG
if (build.model == "yes") {
	source("macs2_model.r")
}
