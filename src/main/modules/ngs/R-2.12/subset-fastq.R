# TOOL subset-fastq.R: "Make a subset of FASTQ" (Returns a subset of N reads from a FASTQ file)
# INPUT reads.fastq: "FASTQ file" TYPE GENERIC
# OUTPUT subset.fastq 
# PARAMETER n.seq: "Size of subset" TYPE INTEGER DEFAULT 100000 (Number of sequences to return from the start of the FASTQ file.)


# binary
binary <- c("head")

# options
n.seq <- as.integer(n.seq)
n.lines <- n.seq * 4

l.opt <- c('-', as.integer(n.lines))
options <- paste(l.opt , collapse = '')

# input files
options <- paste(options,"reads.fastq")

# command
command <- paste(binary, options, " > subset.fastq")

# run
#stop(paste('CHIPSTER-NOTE: ', command))
system(command)
