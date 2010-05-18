# TOOL "Statistics" / find-peaks-using-macs.R: "Find ChIP-seq peaks using MACS" (This tool will search for statistically significantly enriched
# genomic regions in sequencing data from a ChIP-seq experiment. The analysis can be performed on one or more treatment
# samples alone, or relative to one or more control samples.)
# INPUT treatment.tsv: "Treatment data file" TYPE GENERIC
# INPUT control.tsv: "Control data file" TYPE GENERIC
# OUTPUT positive_peaks.tsv: "True enriched peaks"
# OUTPUT analysis_summary.tsv: "Summary of analysis settings and results"
# OUTPUT peak_model.pdf: "A plot of the fitted peak model"
# OUTPUT OPTIONAL negative_peaks.tsv: "The false enriched peaks"
# PARAMETER file.format: "The format of the sequence files" TYPE STRING [ELAND, SAM, BAM, BED] DEFAULT ELAND (The format of the input files.)
# PARAMETER produce.wiggle: "Should wiggle files be produced" TYPE STRING [yes, no] DEFAULT no (Determines if WIGGLE type files should be output or not. By default this option is turned off due to the significantly longer run times it causes. However, for displaying p-values in one track of the Genome Browser, this paramter needs to be yes.)
# PARAMETER species: "The species of the analyzed samples" TYPE STRING [human] DEFAULT human (the species of the samples.)
# PARAMETER read.length: "The length in nucleotides of the sequence reads" TYPE INTEGER FROM 1 TO 200 DEFAULT 30
# PARAMETER band.with: "The scanning window size, typically half the average DNA fragment length" TYPE INTEGER FROM 1 TO 1000 DEFAULT 200
# PARAMETER p.value.threshold: "The unadjusted p-value cutoff for statistical significance" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05
# PARAMETER m.fold: "Sets the m-fold threshold for model building" TYPE INTEGER FROM 1 TO 100 DEFAULT 32

# ALTERNATIVE SETUP
# INPUT sequence{...}.txt: "Sequence data files" TYPE GENERIC
# INPUT phenodata.tsv: "The file containing pehnotype data" TYPE GENERIC
# PARAMETER groups.column: "Column with group labels" TYPE METACOLUMN_SEL DEFAULT groups (Phenodata column describing the experiment groups of the samples. Use "2" for treatment and "1" for control.)
# PARAMETER treatment.group: "The group label used for the treatment samples" TYPE STRING DEFAULT "empty"
# PARAMETER control.group: "The group label used for the control samples" TYPE STRING DEFAULT "empty"

groups.column <- "group"
treatment.group <- "1"
control.group <- "2"
file.format <- "ELAND"
produce.wiggle <- "no"
read.length <- 25
band.with <- 250
p.value.threshold <- 0.05
m.fold <- 32


# find_peals_using_MACS.R
# MG, 17.5.2010

# Loads the libraries!

# Reading data
# files<-dir()
# files<-files[files!="phenodata.tsv"]
# dat<- read.maimages(files=files, columns=columns, annotation=annotation, other.columns=columns.other) 


# Loads the  phenodata files
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
groups<-phenodata[,pmatch(groups.column,colnames(phenodata))]
indices <- seq(1,length(groups),step=1)
treatment.indices <- indices[groups==treatment.group]
if (control.group != "empty") {
	control.indices <- indices[groups==control.group]
}


# Sanity checks
# if(length(unique(groups))==1 | length(unique(groups))>=3) {
#	stop("You need to have exactly two groups to run this analysis")
#}


# If multiple samples per experiment group merge into one single file
if (length(treatment.indices) > 1) {
	command <- "cat"
	command <- paste (command, paste("sequence",treatment.indices[1], sep=""))
	command <- paste(command, ".txt", sep="")
	for (count in 2:length(treatment.indices)) {
		command <- paste(command, " ", sep="")
		command_2 <- paste("sequence", treatment.indices[count], sep="")
		command_2 <- paste(command_2, ".txt", sep="")
		command <- paste(command, command_2, sep="")
	#	assign (paste("data_", count, sep=""), read.table(files[count], header=T, sep="\t")) 
	}
	command <- paste (command, "> treatment.txt")
	system (command)
}
if (length(treatment.indices) == 1) {
	command <- "mv"
	command <- paste(command, paste("sequence", treatment.indices[1], sep=""))
	command <- paste (command, "> treatment.txt")
	system*command()
}
if (length(control.indices) > 1) {
	command <- "cat"
	command <- paste (command, paste("sequence",control.indices[1], sep=""))
	command <- paste(command, ".txt", sep="")
	for (count in 2:length(control.indices)) {
		command <- paste(command, " ", sep="")
		command_2 <- paste("sequence", control.indices[count], sep="")
		command_2 <- paste(command_2, ".txt", sep="")
		command <- paste(command, command_2, sep="")
		#	assign (paste("data_", count, sep=""), read.table(files[count], header=T, sep="\t")) 
	}
	command <- paste (command, "> control.txt")
	system (command)
}
if (length(control.indices) == 1) {
	command <- "mv"
	command <- paste(command, paste("sequence", control.indices[1], sep=""))
	command <- paste (command, "> control.txt")
	system(command)
}

# Trim off any unwanted reads (ambiguous and low quality)
# if (file.format == "ELAND") {
# }

# Remove unmappable reads belonging to random chromosomes or hapmap
system("grep -v random treatment.txt > treatment_2.txt")
system("grep -v hap treatment.txt > treatment_3.txt")
system ("rm -f treatment.txt")
system("rm -f treatment_2.txt")
if (control.group != "empty") {
	system("grep -v random control.txt > control_2.txt")
	system("grep -v hap control.txt > control_3.txt")
	system ("rm -f control.txt")
	system("rm -f control_2.txt")
}

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
		command <- paste("/v/users/chipster/tools/bin/macs", paste(flags, params, collapse=" "))
		if (!is.null(switchOnParams)) {
			switchOnParams <-  paste(switchOnParams, collapse=" ")
			command <- paste(command, switchOnParams)
		}
		# Environment
		environment <- "export PYTHONPATH=/v/users/chipster/tools/lib/python2.6/site-packages ; export PATH=${PATH}:/v/users/akallio/bin ;"
		
#	system ("setenv PYTHONPATH=/v/users/chipster/tools/lib/python2.6/site-packages")
#	setenv PATH=${PATH}:/v/users/chipster/tools/bin
		
		
		# Run macs. Macs writes its output to stderr (stream number 2)
		# &> redirects both stderr and stdout
		# Iterates through mfold values to find low enough that works
		for (mfold in list(32, 24, 16, 8)) {
			system.output <- system(paste(environment, command, paste("--mfold=", mfold, sep=""), "2>", logFile))
#			system.output <- system(paste(command, paste("--mfold=", mfold, sep=""), "2>", logFile))
			if (system.output == 0) {
				break; # was succesfull, don't lower mfold value any more
			}
		}		
		return(invisible(system.output))
	}
}

# Run MACS with default parameters for the data set
runMACS(treatment="./treatment_3.txt", 
		control="control_3.txt", 
		name="./results", 
		format = "ELAND", 
		verbose=3, logFile="./results.log", 
		nomodel=TRUE, help=F, version=FALSE)

# Read in and parse the results

## If final == FALSE, parse the bootsrap results. Else parse the final results
## This feature is not yet implemented.
## <name> is the experiment name and it must be same than in the runMACS
## 10xlog10(pvalue) values are sorted decreasingly"

parseMACSResultsPOS <- function(name, final=FALSE){
	if( final ){
		output <- read.table(file=paste(name, "_peaks.xls", sep=""), skip=0, header=TRUE, stringsAsFactors=FALSE)
		## Choose columns
		output <- output[c("chr","start","end","fold_enrichment", "X.10.log10.pvalue.","tags","summit")]
		## Fix colnames
		colnames(output)[5] <- "-10xlog10pvalue"
		## Sort the result according to the -10xlog10(pvalue)
		output <- output[ order(output[,5], decreasing=TRUE), ]
		return(output)
	}
}


parseMACSResultsNEG <- function(name, final=FALSE){
	if( final ){
		output <- read.table(file=paste(name, "_negative_peaks.xls", sep=""), skip=0, header=TRUE, stringsAsFactors=FALSE)
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
results_TRUE <- parseMACSResultsPOS (name="results",final=TRUE)

## Read in the results for the FALSE, or NEGATIVE, peaks
results_FALSE <- parseMACSResultsNEG (name="results", final=TRUE)


# Write the results to tables to be read into Chipster
write.table(results_TRUE, file="positive-peaks.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(results_FALSE, file="negative-peaks.tsv", sep="\t", quote=FALSE, row.names=FALSE)



