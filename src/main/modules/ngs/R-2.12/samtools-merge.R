# TOOL samtools-merge.R: "Merge BAM files" (Merge BAM files.)
# INPUT alignment{...}.bam: alignment{...}.bam TYPE GENERIC 
# OUTPUT merged.bam 
# OUTPUT merged.bam.bai

# EK 17.6.2011

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# convert sam to bam
system(paste(samtools.binary, "merge merged-not-sorted.bam alignment*.bam"))

# sort bam
system(paste(samtools.binary, "sort merged-not-sorted.bam merged"))

# index bam
system(paste(samtools.binary, "index merged.bam"))

