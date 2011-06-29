# TOOL fastx-clipper.R: "Filter reads for adapters, length and Ns" (Keeps only reads that contain a user-defined adapter sequence. Clips away the adaptor,
# and filters out reads that are too short or contain unknown nucleotides. Adapter-only sequences are removed in the process. This tool is based on the FASTA/Q Clipper tool of the FASTX package. Run on cloud.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT clipped.fastq
# OUTPUT clipped.log
# PARAMETER adapter: "Adapter to be removed" TYPE STRING DEFAULT CCTTAAGG (Adapter sequence that is used for filtering and that is subsequently removed.)
# PARAMETER short: "Discard sequences shorter than" TYPE INTEGER FROM 1 TO 100 DEFAULT 15 (Minimum length of sequences to keep.)



# EK 27.6.2011

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_clipper"))

# command
command <- paste(binary, "-l", short, "-a", adapter, "-i reads.fastq -o clipped.fastq  > clipped.log")

# run
system(command)


