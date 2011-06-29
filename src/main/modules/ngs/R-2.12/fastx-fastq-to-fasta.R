# TOOL fastx-fastq-to-fasta.R: "Convert fastq to fasta" (Convert fastq files to fasta format. This tool is based on the FASTX package.)
# INPUT reads.fastq TYPE GENERIC 
# OUTPUT reads.fasta
# OUTPUT fasta.log
# PARAMETER OPTIONAL remove.unknowns: "Remove sequences with unknown nucleotides" TYPE [yes, no] DEFAULT yes (Remove sequences with unknown nucleotides)
# PARAMETER OPTIONAL rename.identifiers: "Rename sequence identifiers as numbers" TYPE [yes, no] DEFAULT no (Rename sequence identifiers as numbers)



# EK 17.6.2011

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastq_to_fasta"))

# parameters
remove.parameter <- ifelse(remove.unknowns == "yes", "", "-n")
rename.parameter <- ifelse(rename.identifiers == "no", "", "-r")

# command
command <- paste(binary, remove.parameter, rename.parameter, "-i reads.fastq -o reads.fasta > fasta.log")

# run
system(command)

