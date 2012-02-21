# TOOL samtools-subset.R: "Make a subset of BAM" (Retrieves alignments for a specified chromosome or a region, taking the mapping quality into account if needed. For retrieving regions within chromosomes, please use the format chr1:1-100. Please note that in addition to BAM file you have to provide an index file for it. You can create the index file using the tool Index BAM. This tool is based on the SAMtools package.)
# INPUT alignment.bam: "BAM file" TYPE GENERIC 
# INPUT alignment.bai: "Index file .bai" TYPE GENERIC
# OUTPUT alignment-subset.bam
# PARAMETER region: "Region to retrieve alignments for" TYPE STRING DEFAULT chr1 (The genomic region for which you would like to retrieve the alignments for.)
# PARAMETER OPTIONAL mapping.quality: "Minimum mapping quality" TYPE INTEGER FROM 0 TO 1000 DEFAULT 0 (Retrieves only alignments which have a mapping quality higher than this. Note that Bowtie doesn't calculate mapping quality, but just inserts 255 to the field 5 of BAM file if the read aligns, and 0 if it doesn't align.)


# EK 26.10.2011

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

system(paste(samtools.binary, "view -b -q", mapping.quality, "-o alignment-subset.bam alignment.bam", region))


#test.command <- paste(samtools.binary, "view -b -q", mapping.quality, "-o alignment-subset.bam alignment.bam", region)
#stop(paste('CHIPSTER-NOTE: ', test.command))




