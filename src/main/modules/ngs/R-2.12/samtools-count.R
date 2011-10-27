# TOOL samtools-count.R: "Count alignments in BAM file" (Counts alignments in a BAM file, taking the mapping quality into account if needed. This tool is based on the SAMtools package.)
# INPUT alignment.bam TYPE GENERIC 
# OUTPUT alignment-counts.txt 
# PARAMETER OPTIONAL quality: "Minimum mapping quality required for counting" TYPE INTEGER FROM 0 TO 1000 DEFAULT 0 (Counts only alignments which have a mapping quality higher than this. Note that Bowtie doesn't calculate mapping quality, but just inserts 255 to the field 5 of BAM file if the read aligns, and 0 if it doesn't align.)


# EK 26.10.2011

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

system(paste(samtools.binary, "view -bc -q", mapping.quality, "-o alignment-counts.txt alignment.bam"))



