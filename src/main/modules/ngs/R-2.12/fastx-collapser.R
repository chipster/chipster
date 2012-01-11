# TOOL fastx-collapser.R: "Remove duplicate reads from FASTQ" (Collapses identical sequences in a FASTQ file into a single sequence. The sequences are renamed with sequence number and the multiplicity value. This tool is based on the FASTA/Q Collapser tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT duplicates-removed.fastq 


# EK 12.1.2012

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_collapser"))

# command
command <- paste(binary, "-i reads.fastq -o duplicates-removed.fastq")

# run
system(command)