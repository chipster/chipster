# TOOL samtools-sort-index-BAM.R: "Sort and index BAM" (Sorts BAM by chromosomal location and builds an index for it. This tool is based on the SAMtools package.)
# INPUT alignment.bam: "BAM file to sort and index" TYPE GENERIC 
# OUTPUT alignment-sorted.bam 
# OUTPUT alignment-sorted.bam.bai 

# EK 10.1.2012

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# sort bam
system(paste(samtools.binary, "sort alignment.bam alignment-sorted"))

# index bam
system(paste(samtools.binary, "index alignment-sorted.bam"))

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("alignment-sorted.bam", paste(inputnames$alignment.bam))
outputnames[2,] <- c("alignment-sorted.bam.bai", paste(inputnames$alignment.bam, ".bai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)