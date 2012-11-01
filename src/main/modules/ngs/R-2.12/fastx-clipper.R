# TOOL fastx-clipper.R: "Filter reads for adapters, length and Ns" (Filters reads accordingto user-defined adapter sequence. Clips away the adapter,
# and optionally filters out reads that are too short after clipping or contain unknown nucleotides. Adapter-only sequences are removed in the process. This tool is based on the FASTA/Q Clipper tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT clipped.fastq
# OUTPUT clipped.log
# PARAMETER adapter: "Adapter to be removed" TYPE STRING DEFAULT CCTTAAGG (Adapter sequence that is used for filtering and that is subsequently removed.)
# PARAMETER minimum.alignment: "Minimum adapter alignment length" TYPE INTEGER FROM 0 DEFAULT 0 (Required minimum adapter alignment length. Maximum is adapter length. 0 means option is ignored.)
# PARAMETER short: "Discard sequences shorter than" TYPE INTEGER FROM 1 DEFAULT 15 (Minimum length of sequences to keep.)
# PARAMETER quality.format: "Quality value format used" TYPE [sanger: Sanger, illuminaold: "Illumina GA v1.3-1.5"] DEFAULT sanger (What quality encoding is used in your FASTQ file. Select Sanger if your data comes from Illumina 1.8 or later, SOLiD or 454.)
# PARAMETER discard.n: "Discard sequences with Ns" TYPE [yes, no] DEFAULT yes (Keep sequences with unknown nucleotides. Default is to discard such sequences.)
# PARAMETER output.options: "Output options" TYPE [clipped: "Keep only clipped reads", unclipped: "Keep only non-clipped reads", both: "Keep both clipped and non-clipped reads"] DEFAULT clipped (You can choose to keep only clipped reads (reads that contained adapter\), only non-clipped reads (reads that did not contain an adapter\) or both clipped and non-clipped reads.)

# EK 27.6.2011
# AMS 01.11.2012: Added minimum.alignement and output.options options

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_clipper"))

# options
options <- paste("")
options <- paste(options, "-a", adapter)
if (minimum.alignment > 0) {
	options <- paste(options, "-M", minimum.alignment)
}
options <- paste(options, "-l", short)
if (quality.format == "sanger") {
	options <- paste(options, "-Q 33")
}
if (discard.n == "no") {
	options <- paste(options, "-n")
}
if (output.options == "clipped") {
	options <- paste(options, "-c")
}
if (output.options == "unclipped") {
	options <- paste(options, "-C")
}

# command
command <- paste(binary, "-v", options, "-i reads.fastq -o clipped.fastq  > clipped.log")

# run
system(command)


