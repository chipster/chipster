# TOOL samtools-idxstats.R: "Statistics for a BAM file" (Counts how many alignments there are per each chromosome. Please note that in addition to BAM file you have to provide an index file with the same name as the BAM. You can create the index file using the tool Convert SAM to BAM, sort and index BAM. This tool is based on the SAMtools package.)
# INPUT alignment.bam: "BAM file" TYPE GENERIC 
# INPUT alignment.bai: "Index file .bai" TYPE GENERIC
# OUTPUT bam-stats.tsv

# EK 26.10.2011

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# convert sam to bam
system(paste(samtools.binary, "idxstats alignment.bam > bam-stats.tsv"))

