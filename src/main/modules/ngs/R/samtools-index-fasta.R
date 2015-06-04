# TOOL samtools-index-fasta.R: "Index FASTA" (Creates an index for a FASTA file. This tool is based on the SAMTools package.)
# INPUT sequence.fa TYPE GENERIC 
# OUTPUT sequence.fa.fai

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequence.fa")


# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

system(paste(samtools.binary, "faidx sequence.fa"))

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("sequence.fa.fai", paste(remove_postfix(inputnames$sequence.fa, ".gz"), ".fai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)

