# TOOL fastx-fastq-to-fasta.R: "Convert FASTQ to FASTA" (Convert FASTQ files to FASTA format. This tool is based on the FASTX package.)
# INPUT reads.fastq TYPE GENERIC 
# OUTPUT reads.fasta
# OUTPUT fasta.log
# PARAMETER OPTIONAL remove.unknowns: "Remove sequences with unknown nucleotides" TYPE [yes, no] DEFAULT yes (Remove sequences with unknown nucleotides)
# PARAMETER OPTIONAL rename.identifiers: "Rename sequence identifiers as numbers" TYPE [yes, no] DEFAULT no (Rename sequence identifiers as numbers)
# PARAMETER quality.format: "Quality value format used" TYPE [sanger: Sanger, illuminaold: "Illumina GA v1.3-1.5"] DEFAULT sanger (What quality encoding is used in your FASTQ file. Select Sanger if your data comes from Illumina 1.8 or later, SOLiD or 454.)



# EK 17.6.2011

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastq_to_fasta"))

# parameters
remove.parameter <- ifelse(remove.unknowns == "yes", "", "-n")
rename.parameter <- ifelse(rename.identifiers == "no", "", "-r")
quality.scale <- ifelse(quality.format == "sanger", "-Q 33", "")

# command
command <- paste(binary, remove.parameter, rename.parameter, quality.scale, "-v -i reads.fastq -o reads.fasta > fasta.log")

# run
system(command)

