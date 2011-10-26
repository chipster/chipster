# TOOL fastx-renamer.R: "Rename reads" (Renames read identifiers in a FASTQ file. The original names can be replaced by a running number or by the sequence itself. This tool is based on the Renamer tool of the FASTX package.) 
# INPUT reads.fastq TYPE GENERIC 
# OUTPUT renamed-reads.fastq
# PARAMETER rename.type: "Rename read identifiers to" TYPE [COUNT: "Running number", SEQ: "Sequence itself"] DEFAULT COUNT (Should the reads be named with running number or the sequence itself.)


# EK 24.10.2011

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_renamer"))

# command
command <- paste(binary, "-n", rename.type, "-i reads.fastq -o renamed-reads.fastq")

# run
system(command)

