# TOOL fastx-collapser.R: "Remove duplicate reads from FASTQ or FASTA" (Collapses identical sequences in a FASTQ/FASTA file into a single sequence. The sequences are renamed with sequence number and the multiplicity value. This tool is based on the FASTA/Q Collapser tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT duplicates-removed.fastq 
# PARAMETER quality.format: "Quality value format used" TYPE [sanger: Sanger, illuminaold: "Illumina GA v1.3-1.5"] DEFAULT sanger (What quality encoding is used in your FASTQ file. Select Sanger if your data comes from Illumina 1.8 or later, SOLiD or 454.)



# EK 12.1.2012

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_collapser"))

# command
quality.scale <- ifelse(quality.format == "sanger", "-Q 33", "")
command <- paste(binary, quality.scale, "-i reads.fastq -o duplicates-removed.fastq")

# run
system(command)