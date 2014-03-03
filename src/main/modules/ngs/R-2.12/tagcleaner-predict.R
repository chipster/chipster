# TOOL tagcleaner-predict.R: "Predict primers/adaptors with TagCleaner" (Tool will attempt to predict the tag (i.e. the primer/adapter\) at either or both sites, if possible. The algorithm implemented for the tag prediction assumes the randomness of a typical metagenome. Datasets that do not contain random sequences from organisms in an environment, but rather contain, for example, 16S data may cause incorrect detection of the tag sequences. However, the tag sequences will most likely be over-predicted and can be redefined by the user prior to data processing. The tag sequence prediction uses filtered base frequencies instead of raw base frequencies. This allows a more accurate prediction as it accounts for incomplete and shifted tag sequences. If no tags are reported, then no tags could be identified in the data set. No trimming will be performed. This tool is based on TagCleaner.)
# INPUT reads: "FASTQ/FASTA file" TYPE GENERIC
# OUTPUT OPTIONAL tag.predict.tsv
# PARAMETER input.type: "Input type" TYPE [FASTQ, FASTA] DEFAULT FASTQ (Input type.)

# AMS 2013.02.18

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads")

# binary
tagcleaner.binary <- c(file.path(chipster.tools.path, "tagcleaner", "tagcleaner.pl"))

# options
# options for input type
options <- paste("")
if (input.type == "FASTQ"){
	options <- paste(options, "-fastq reads")
}
if (input.type == "FASTA"){
	options <- paste(options, "-fasta reads")
}
# common options
options <- paste(options, "-64")
options <- paste(options, "-predict")

# command
command <- paste(tagcleaner.binary, options, "> tag.predict.tsv")

#stop(paste('CHIPSTER-NOTE: ', command))
# run
system(command)