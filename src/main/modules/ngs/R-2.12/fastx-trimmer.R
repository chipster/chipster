# TOOL fastx-trimmer.R: "Trim reads" (Trims reads to a user-specified length. This tool is based on the FASTA/Q Trimmer tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT trimmed.fastq 
# PARAMETER first: "First base to keep" TYPE INTEGER FROM 1 TO 100 DEFAULT 1 (First base to keep.)
# PARAMETER last: "Last base to keep" TYPE INTEGER FROM 1 TO 300 DEFAULT 75 (Last base to keep.)




# EK 17.6.2011

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_trimmer"))

# command
command <- paste(binary, "-f", first, "-l", last, "-i reads.fastq -o trimmed.fastq")

# run
system(command)