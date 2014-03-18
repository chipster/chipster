# TOOL samtools-index.R: "Index BAM" (Creates an index for a BAM file. Please rename the index file to have the same name as the BAM file. This tool is based on the SAMtools package.)
# INPUT alignment.bam TYPE GENERIC 
# OUTPUT alignment.bam.bai


# EK 30.10.2011

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

system(paste(samtools.binary, "index alignment.bam > alignment.bam.bai"))






