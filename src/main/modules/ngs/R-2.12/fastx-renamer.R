# TOOL fastx-renamer.R: "Rename reads" (Renames read identifiers in a FASTQ file. The original names can be replaced by a running number or by the sequence itself. This tool is based on the Renamer tool of the FASTX package.) 
# INPUT reads.fastq TYPE GENERIC 
# OUTPUT renamed-reads.fastq
# PARAMETER rename.type: "Rename read identifiers to" TYPE [COUNT: "Running number", SEQ: "Sequence itself"] DEFAULT COUNT (Should the reads be named with running number or the sequence itself.)
# PARAMETER quality.format: "Quality value format used" TYPE [sanger: Sanger, illuminaold: "Illumina GA v1.3-1.5"] DEFAULT sanger (What quality encoding is used in your FASTQ file. Select Sanger if your data comes from Illumina 1.8 or later, SOLiD or 454.)


# EK 24.10.2011

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_renamer"))

# command
quality.scale <- ifelse(quality.format == "sanger", "-Q 33", "")
command <- paste(binary, "-n", rename.type, quality.scale, "-i reads.fastq -o renamed-reads.fastq")

# run
system(command)

