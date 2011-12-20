# TOOL samtools-idxstats.R: "Statistics for BAM" (Counts how many alignments there are per each chromosome. Please note that in addition to BAM file you have to provide an index file with the same name. You can create the index file using the tool Index BAM. This tool is based on the SAMtools package.)
# INPUT alignment.bam: "BAM file" TYPE GENERIC 
# INPUT alignment.bai: "Index file .bai" TYPE GENERIC
# OUTPUT bam-stats.tsv

# EK 26.10.2011

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# convert sam to bam
system(paste(samtools.binary, "idxstats alignment.bam > bam-out.tsv"))

# read into R to prettify table
dat <- read.table(file="bam-out.tsv", header=FALSE, sep="\t")
colnames (dat) [1] <- "chr"
colnames (dat) [2] <- "chr_length"
colnames (dat) [3] <- "reads_mapped"
colnames (dat) [4] <- "reads_unmapped"

# output table with headers
write.table(dat, file="bam-stats.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# EOF
