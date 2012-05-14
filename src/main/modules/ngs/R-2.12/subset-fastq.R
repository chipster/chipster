# TOOL subset-fastq.R: "Subset FASTQ" (Returns a subset of N sequences from a FASTQ file)
# INPUT reads.fastq: "FASTQ file" TYPE GENERIC
# OUTPUT subset.fastq 
# PARAMETER n.seq: "Size of subset" TYPE INTEGER DEFAULT 100000 (Number of sequences to return from the start of the FASTQ file.)


# binary
binary <- c("head")

# options
n.seq <- as.integer(n.seq)
n.lines <- n.seq * 4

options <- paste(' ', '-', n.lines, collapse = '')

# input files
options <- paste(options,"reads.fastq")

# command
command <- paste(binary, options, " > subset.fastq")

# run
system(command)
