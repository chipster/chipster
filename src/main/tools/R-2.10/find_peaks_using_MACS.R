# ANALYSIS Statistics/"Find ChIP-seq peaks using MACS" (This tool will search for statistically significantly enriched
# genomic regions in sequencing data from a ChIP-seq experiment. The analysis can be performed on one or more treatment
# samples alone, or relative to one or more control samples.)
# INPUT SHORT_READS sequence[...].tsv
# OUTPUT positive_peaks.tsv, negative_peaks.tsv, analysis_summary.tsv, peak_model.pdf
# PARAMETER groups.column METACOLUMN_SEL DEFAULT group (Phenodata column describing the experiment groups of the samples. Use "2" for treatment and "1" for control.)
# PARAMETER file.format [ELAND, BED, SAM, BAM] DEFAULT ELAND (The format of the sequence files.)
# PARAMETER produce.WIGGLE [yes, no] DEFAULT no (Determines if WIGGLE type files should be output or not.
# By default this option is turned off due to the significantly longer run times it causes. However, for
# displaying p-values in one track of the Genome Browser, this paramter needs to be "yes".)
# PARAMETER read.length
# PARAMETER band.with
# PARAMETER p.value.threshold
# PARAMETER m.fold

# find_peals_using_MACS.R
# MG, 22.4.2010

# Loads the libraries


# Loads the normalized data and phenodata files
data_1 <- read.table(file="normalized_mirna.tsv", header=T, sep="\t", row.names=1)
phenodata_1 <- read.table("phenodata_mirna.tsv", header=T, sep="\t")

# If multiple samples per experiment group merge into one single file

# Trim off any unwanted reads (ambiguous and low quality)

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
		command <- paste("macs", paste(flags, params, collapse=" "))
		if (!is.null(switchOnParams)) {
			switchOnParams <-  paste(switchOnParams, collapse=" ")
			command <- paste(command, switchOnParams)
		}
		# Environment
		environment <- "export PYTHONPATH /v/users/chipster/tools/lib/python2.6/site-packages ; export PATH ${PATH}:/v/users/chipster/tools/bin ; "
		# Run macs. Macs writes its output to stderr (stream number 2)
		# &> redirects both stderr and stdout
		# Iterates through mfold values to find low enough that works
		for (mfold in list(32, 24, 16, 8)) {
#			system.output <- system(paste(environment, command, paste("--mfold=", mfold, sep=""), "2>", logFile))
			system.output <- system(paste(command, paste("--mfold=", mfold, sep=""), "2>", logFile))
			if (system.output == 0) {
				break; # was succesfull, don't lower mfold value any more
			}
		}		
		return(invisible(system.output))
	}
}

# Run MACS with default parameters for FoxA1 test data set
runMACS(treatment="/fs/local/users/chipster/tools/Seq_data/FoxA1/Input_tags.bed", 
		control="/fs/local/users/chipster/tools/Seq_data/FoxA1/Treatment_tags.bed", 
		name="/fs/local/users/chipster/tools/Seq_data/FoxA1/MACS_results/FoxA1_t30_b175_p1e5_m10", 
		format = "BED", 
		verbose=3, logFile="/fs/local/users/chipster/tools/Seq_data/FoxA1/MACS_results/FoxA1_t30_b175_p1e5_m10.log", 
		nomodel=TRUE, help=F, version=FALSE)

# Read in and parse the results

## If final == FALSE, parse the bootsrap results. Else parse the final results
## This feature is not yet implemented.
## <name> is the experiment name and it must be same than in the runMACS
## 10xlog10(pvalue) values are sorted decreasingly"

parseMACSResultsPOS <- function(name, final=FALSE){
	if( final ){
		output <- read.table(file=paste(name, "_peaks.xls", sep=""), skip=3, header=TRUE, stringsAsFactors=FALSE)
		## Choose columns
		output <- output[c("chr","start","end","fold_enrichment", "X.10.log10.pvalue.","tags","summit")]
		## Fix colnames
		colnames(output)[5] <- "-10xlog10pvalue"
		## Sort the result according to the -10xlog10(pvalue)
		output <- output[ order(output[,5], decreasing=TRUE), ]
		return(output)
	}
}

## Read in the results for the TRUE peaks
results_TRUE <- parseMACSResultsPOS (name="/fs/local/users/chipster/tools/Seq_data/FoxA1/MACS_results/FoxA1_t30_b175_p1e5_m10",
		final=TRUE)

## Read in the results for the FALSE, or NEGATIVE, peaks
results_TRUE <- parseMACSResultsPOS (name="/fs/local/users/chipster/tools/Seq_data/FoxA1/MACS_results/FoxA1_t30_b175_p1e5_m10",
		final=TRUE)

