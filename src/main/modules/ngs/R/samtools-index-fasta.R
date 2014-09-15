# TOOL samtools-index-fasta.R: "Index FASTA" (Creates an index for a FASTA file. Please rename the index file to have the same name as the FASTA file. This tool is based on the SAMTools package.)
# INPUT sequence.fa TYPE GENERIC 
# OUTPUT sequence.fa.fai

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

system(paste(samtools.binary, "faidx sequence.fa"))






