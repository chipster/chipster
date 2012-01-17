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
