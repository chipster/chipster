# TOOL dexseq-exoncount.R: "Count aligned reads per exons for DEXSeq" (Given mapped reads in a sorted BAM file, this tool counts the reads that fall into each non-overlapping exonic part using the script dexseq-count.py. In order to use the output in DEXSeq, you need to select all samples and run the tool \"Utilities - Define NGS experiment\".)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# OUTPUT exon-counts.tsv
# PARAMETER paired: "Does the alignment file contain paired-end data" TYPE [yes, no] DEFAULT no (Does the alignment data contain paired end or single end reads?)
# PARAMETER sort: "Sorting order of the BAM file" TYPE [pos: coordinate, name: name] DEFAULT pos (Were reads sorted by read names or chromosomal coordinates? Relevant for paired end data only.)
# PARAMETER organism: "Organism" TYPE [Homo_sapiens.GRCh37.68.chr.DEXSeq.gtf: "Human (hg19.68)", Mus_musculus.GRCm38.68.chr.DEXSeq.gtf: "Mouse (mm10.68)", Rattus_norvegicus.RGSC3.4.68.chr.DEXSeq.gtf: "Rat (rn4.68)"] DEFAULT Homo_sapiens.GRCh37.68.chr.DEXSeq.gtf (Which organism is your data from.)

# 18.09.2012 TH and EK 
# 16.07.2013 EK, BAM sorting changed
# 23.04.2013 MK, polished script and added parameter defining the sort order of paired end data

# convert bam to sam, sort bam if the data is paired-end
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
samtools.view <- paste(samtools.binary, "view alignment.bam")

# gtf
gtf <- file.path(chipster.tools.path, "genomes", "gtf", organism)

# exoncount
dexseq.binary <- file.path(chipster.tools.path, "dexseq-exoncounts", "dexseq_count.py")
if(paired == "yes") {
	paired.end.data <- ifelse(sort == "pos", paste("-p yes -r pos"), paste("-p yes -r name"))
}
dexseq.command <- paste("python", dexseq.binary, paired.end.data, gtf,"- exon-counts.tsv")

# run 
command <- paste(samtools.view, "|", dexseq.command)
system(command)

# bring in file to R environment for formating
file <- c("exon-counts.tsv")
dat <- read.table(file, header=F, sep="\t")
names(dat) <- c("id", "count")

# write result table to output
write.table(dat, file="exon-counts.tsv", col.names=T, quote=F, sep="\t", row.names=F)


# EOF


