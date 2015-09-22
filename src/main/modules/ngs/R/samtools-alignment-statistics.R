# TOOL samtools-alignment-statistics.R: "Alignment statistics with Samtools" (Counts statistics for alignments in a BAM file. This tool is based on the SAMtools package.)
# INPUT alignment.bam TYPE GENERIC 
# OUTPUT alignment-statistics.txt 

# AMS 17.09.2015

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

system(paste(samtools.binary, "flagstat", "alignment.bam > alignment-statistics.txt"))






