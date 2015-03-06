# TOOL dexseq-exoncount-own-gtf.R: "Count aligned reads per exons for DEXSeq using own GTF" (Given mapped reads in a BAM file, this tool counts the reads that fall into each non-overlapping exonic part using the script dexseq-count.py. In order to use the output in DEXSeq, you need to select all samples and run the tool \"Utilities - Define NGS experiment\".)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# INPUT annotation.gtf: "GTF annotation file" TYPE GENERIC
# OUTPUT exon-counts.tsv
# OUTPUT OPTIONAL exon-counts-info.txt
# PARAMETER paired: "Does the BAM file contain paired-end data" TYPE [yes, no] DEFAULT no (Does the alignment data contain paired end or single end reads?)
# PARAMETER stranded: "Was the data produced with a strand-specific protocol" TYPE [yes, no, reverse] DEFAULT no (Select no if your data was not produced with a strand-specific RNA-seq protocol, so that a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature. If you select yes, the read has to be mapped to the same strand as the feature.)
# PARAMETER OPTIONAL mode: "Mode to handle reads overlapping more than one feature" TYPE [union, intersection-strict, intersection-nonempty] DEFAULT union (How to deal with reads that overlap more than one gene or exon?)

# 18.09.2012 TH and EK 
# 16.07.2013 EK, BAM sorting changed
# 23.04.2013 MK, added the info output file and strandedness parameter
# 01.06.2014 EK, fixed BAM sorting by name, updated to use dexseq-count.py from DEXSeq v1.8.0, added support for BAMs which don't have the chr prefix in chromosome names, moved NH tag production to a separate script
# 03.06.2014 AMS, changed the way chr is handled, updated gtf files

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("annotation.gtf")

# if BAM contains paired-end data, sort it by read names
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
if(paired == "yes"){
	system(paste(samtools.binary, "sort -n alignment.bam name-sorted"))
	bam<-"name-sorted.bam"
} else {
	bam<-"alignment.bam"
}

# User provided GTF need to be preprocessed
prepare.binary <- file.path(chipster.tools.path, "dexseq-exoncounts", "dexseq_prepare_annotation.py")
prepare.command <- paste("python", prepare.binary, "annotation.gtf annotation.dexseq.gtf")
system(prepare.command)

# counts reads per non-overlapping exonic regions
dexseq.binary <- file.path(chipster.tools.path, "dexseq-exoncounts", "dexseq_count.py")
dexseq.command <- paste("python", dexseq.binary, "-f bam -r name -s", stranded, "-p", paired, "annotation.dexseq.gtf", bam, "exon-counts-out.tsv")
system(dexseq.command)

# separate result file
system("head -n -4 exon-counts-out.tsv > exon-counts.tsv")
system("tail -n 4 exon-counts-out.tsv > exon-counts-info.txt")

# bring in file to R environment for formating
file <- c("exon-counts.tsv")
dat <- read.table(file, header=F, sep="\t")
names(dat) <- c("id", "count")

# write result table to output
write.table(dat, file="exon-counts.tsv", col.names=T, quote=F, sep="\t", row.names=F)

# EOF

