# TOOL samtools-rmdup.R: "Remove duplicate reads from BAM" (Removes potential PCR duplicates. You have to indicate if your data is single end or paired end. In the paired end mode the ISIZE has to be correctly set. This tool is based on the SAMtools package.)
# INPUT alignment.bam TYPE GENERIC 
# OUTPUT duplicates-removed.bam
# PARAMETER end.type: "Is the data paired end or single end reads" TYPE [paired, single] DEFAULT paired (Does the data come from paired end or single end sequencing?)


# EK 12.1.2012

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# parameters
single.end <- ifelse(end.type == "paired", "", "-s")

# command
system(paste(samtools.binary, "rmdup", single.end, "alignment.bam duplicates-removed.bam"))

