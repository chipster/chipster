# TOOL samtools-index.R: "Index BAM" (Creates an index for a BAM file. This tool is based on the SAMtools package.)
# INPUT alignment.bam TYPE GENERIC 
# OUTPUT alignment.bam.bai


# EK 30.10.2011

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

system(paste(samtools.binary, "index alignment.bam > alignment.bam.bai"))

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("alignment.bam.bai", paste(inputnames$alignment.bam, ".bai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)






