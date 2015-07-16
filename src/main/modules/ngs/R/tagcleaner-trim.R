# TOOL tagcleaner-trim.R: "Trim primers/adaptors with TagCleaner" (Trim tags (i.e. the primers/adapters\) from sequences. This tool is based on TagCleaner.)
# INPUT reads: "FASTQ/FASTA file" TYPE GENERIC
# INPUT OPTIONAL qual: "QUAL file" TYPE GENERIC
# OUTPUT trim.log.txt
# OUTPUT OPTIONAL trimmed.fastq.gz
# OUTPUT OPTIONAL trimmed.fasta.gz
# OUTPUT OPTIONAL trimmed.qual.gz
# PARAMETER input.type: "Input type" TYPE [FASTQ, FASTA] DEFAULT FASTQ (Input type. If using a FASTA file you can provide an optional QUAL file to trim sequences and their quality scores.)
# PARAMETER tag5: "Tag sequence at 5'-end" TYPE STRING DEFAULT "-" (Tag sequence at 5'-end.)
# PARAMETER tag3: "Tag sequence at 3'-end" TYPE STRING DEFAULT "-" (Tag sequence at 3'-end.)
# PARAMETER OPTIONAL mm5: "Maximum  mismatches at the 5'-end" TYPE INTEGER DEFAULT 0 (Maximum number of allowed mismatches at the 5'-end.)
# PARAMETER OPTIONAL mm3: "Maximum  mismatches at the 3'-end" TYPE INTEGER DEFAULT 0 (Maximum number of allowed mismatches at the 3'-end.)
# PARAMETER OPTIONAL trimwithin: "Trim within" TYPE INTEGER DEFAULT 0 (The sequence of the tag could occur not only at the sequence end, but also at any other position of the sequence. To assure that only tags are trimmed, the tag sequences can be defined to occur only at the ends allowing a certain number of variable bases. If specified, the value has to be at least the number of bases in the tag sequence. 0 value means option is ignored.)
# PARAMETER OPTIONAL split: "Split fragment-to-fragment concatenations" TYPE [yes, no] DEFAULT yes (This option removes tag contaminations inside the sequences and splist fragment-to-fragment concatenations into separate sequences. This feature should be used with caution for inputs with only a 5' or 3' tag sequence (likely splits too many false positive that naturally occur for single tags compared to much longer concatenated 5' and 3' tags\).)
# PARAMETER OPTIONAL split.mismatch: "Allowed mismatches" TYPE INTEGER DEFAULT 0 (Maximum number of allowed mismatches for the internal (concatenated\) tag sequence(s\).)

# AMS 2013.02.18
# AMS 11.3.2014, gzip fastq outputs

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads")
unzipIfGZipFile("qual")

# binary
tagcleaner.binary <- c(file.path(chipster.tools.path, "tagcleaner", "tagcleaner.pl"))

# options
options <- paste("")
# options for input type
if (input.type == "FASTQ"){
	options <- paste(options, "-fastq reads")
}
if (input.type == "FASTA"){
	options <- paste(options, "-fasta reads")
	if (file.exists("qual")){
		options <- paste(options, "-qual qual")
	}
}
# options for parameters
notags <- TRUE
if (tag3 != "-"){
	options <- paste(options, "-tag3", tag3)
	options <- paste(options, "-mm3", mm3)
	notags <- FALSE	
}
if (tag5 != "-"){
	options <- paste(options, "-tag5", tag5)
	options <- paste(options, "-mm5", mm5)
	notags <- FALSE
}
if (notags){
	stop('CHIPSTER-NOTE: No tag sequence specified. If tags are unknow, you can first use the "Predict primers/adaptors" tool.')
}
if (trimwithin > 0){
	if (trimwithin < max(nchar(tag3),nchar(tag5))){
		stop('CHIPSTER-NOTE: Value for "Allowed variable bases" has to be at least the number of bases in the tag sequence.')
	}
	options <- paste(options, "-trim_within", trimwithin)
}
if (split == "yes"){
	options <- paste(options, "-split", split.mismatch)
}
# common options
options <- paste(options, "-64")
options <- paste(options, "-verbose -log trim.log.txt")
options <- paste(options, "-out trimmed")

# command
command <- paste(tagcleaner.binary, options)

# run
system(command)
system("gzip *.fastq")
system("gzip *.fasta")
system("gzip *.qual")

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow=3, ncol=2)

base1 <- strip_name(inputnames$reads)

outputnames[1,] <- c("trimmed.fastq.gz", paste(base1, ".fastq.gz", sep =""))
outputnames[2,] <- c("trimmed.fasta.gz", paste(base1, ".fasta.gz", sep =""))
outputnames[3,] <- c("trimmed.qual.gz", paste(base1, ".qual.gz", sep =""))


# Write output definitions file
write_output_definitions(outputnames)
