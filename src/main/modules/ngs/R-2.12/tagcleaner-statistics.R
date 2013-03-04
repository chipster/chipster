# TOOL tagcleaner-statistics.R: "Statistics for primers/adaptors" (Prints the number of tag sequences matching for different numbers of mismatches. No trimming will be performed. This tool is based on TagCleaner.)
# INPUT reads: "FASTQ/FASTA file" TYPE GENERIC
# INPUT OPTIONAL qual: "QUAL file" TYPE GENERIC
# OUTPUT tag.statistics.tsv
# PARAMETER input.type: "Input type" TYPE [FASTQ, FASTA+QUAL] DEFAULT FASTQ (Input type. If using a FASTA file you can also optionally provide a QUAL file.)
# PARAMETER tag3: "Tag sequence at 3'-end" TYPE STRING DEFAULT "-" (Tag sequence at 3'-end.)
# PARAMETER tag5: "Tag sequence at 5'-end" TYPE STRING DEFAULT "-" (Tag sequence at 5'-end.)


# AMS 2013.02.18

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads")

# binary
tagcleaner.binary <- c(file.path(chipster.tools.path, "tagcleaner", "tagcleaner.pl"))

# options
options <- paste("")
if (input.type == "FASTQ"){
	options <- paste(options, "-fastq reads")
}
if (input.type == "FASTA+QUAL"){
	options <- paste(options, "-fasta reads")
	if (file.exists("qual")){
		options <- paste(options, "-qual qual")
	}
}
# options for parameters
notags <- TRUE
if (tag3 != "-"){
	options <- paste(options, "-tag3", tag3)
	notags <- FALSE	
}
if (tag5 != "-"){
	options <- paste(options, "-tag5", tag5)
	notags <- FALSE
}
if (notags){
	stop('CHIPSTER-NOTE: No tag sequence specified. If tags are unknow, you can first use the "Predict primers/adaptors" tool.')
}
# common options
options <- paste(options, "-64")
options <- paste(options, "-stats")

# command
command <- paste(tagcleaner.binary, options, "> tag.statistics.tsv")

# run
system(command)