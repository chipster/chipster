# TOOL fastx-trimmer.R: "Trim reads" (Trims reads to a user-specified length. This tool is based on the FASTA/Q Trimmer tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT trimmed.fastq 
# PARAMETER first: "First base to keep" TYPE INTEGER FROM 1 TO 100 DEFAULT 1 (First base to keep.)
# PARAMETER last: "Last base to keep" TYPE INTEGER FROM 1 TO 300 DEFAULT 75 (Last base to keep.)
# PARAMETER quality.format: "Quality value format used" TYPE [sanger: Sanger, illuminaold: "Illumina GA v1.3-1.5"] DEFAULT sanger (What quality encoding is used in your FASTQ file. Select Sanger if your data comes from Illumina 1.8 or later, SOLiD or 454.)




# EK 17.6.2011

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_trimmer"))

# command
quality.scale <- ifelse(quality.format == "sanger", "-Q 33", "")
command <- paste(binary, "-f", first, "-l", last, quality.scale, "-i reads.fastq -o trimmed.fastq")

# run
system(command)