# TOOL samtools-convert-SAMtoBAM-sort-index.R: "Convert SAM to BAM, sort and index" (Converts SAM file to BAM, sorts it by chromosomal location, and builds an index for it. This tool is based on the SAMtools package.)
# INPUT alignment.sam: "SAM file to convert" TYPE GENERIC 
# OUTPUT alignment.bam 
# OUTPUT alignment.bam.bai 

# EK 10.1.2012

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# convert sam to bam
system(paste(samtools.binary, "view -bS alignment.sam -o alignmentUnsorted.bam"))

# sort bam
system(paste(samtools.binary, "sort alignmentUnsorted.bam alignment"))

# index bam
system(paste(samtools.binary, "index alignment.bam"))

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

base <- strip_name(inputnames$alignment.sam)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("alignment.bam", paste(base, ".bam", sep =""))
outputnames[2,] <- c("alignment.bam.bai", paste(base, ".bam.bai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)