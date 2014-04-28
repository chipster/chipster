# TOOL dexseq-exoncount.R: "Count aligned reads per exons for DEXSeq" (Given mapped reads in a sorted BAM file, this tool counts the reads that fall into each non-overlapping exonic part using the script dexseq-count.py. In order to use the output in DEXSeq, you need to select all samples and run the tool \"Utilities - Define NGS experiment\".)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# OUTPUT exon-counts.tsv
# OUTPUT OPTIONAL exon-counts-info.tsv
# PARAMETER paired: "Does the alignment file contain paired-end data" TYPE [yes, no] DEFAULT no (Does the alignment data contain paired end or single end reads?)
# PARAMETER sort: "Sorting order of the BAM file" TYPE [pos: coordinate, name: name] DEFAULT pos (Were reads sorted by read names or chromosomal coordinates? Relevant for paired end data only.)
# PARAMETER organism: "Organism" TYPE [Homo_sapiens.GRCh37.68.chr.DEXSeq.gtf: "Human (hg19.68)", Mus_musculus.GRCm38.68.chr.DEXSeq.gtf: "Mouse (mm10.68)", Rattus_norvegicus.RGSC3.4.68.chr.DEXSeq.gtf: "Rat (rn4.68)"] DEFAULT Homo_sapiens.GRCh37.68.chr.DEXSeq.gtf (Which organism is your data from.)
# PARAMETER OPTIONAL stranded: "Was the data produced with a strand-specific RNA-seq protocol" TYPE [yes, no, reverse] DEFAULT no (If you select no, a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature. If you select yes, the read has to be mapped to the same strand as the feature. You have to say no, if your was not made with a strand-specific RNA-seq protocol, because otherwise half your reads will be lost.)
# PARAMETER OPTIONAL addnh: "Add NH flag" TYPE [yes: yes, no: no] DEFAULT no (Extract uniquely matched reads through MAPQ field, i.e. reads with MAPQ higher than 3, and append alignment lines by NH:i:1 flags. Please note that MAPQ works only if MAPQ values are properly set and that filtering also removes low-quality reads. Preferably, this kind of filtering should be carried out already at the mapping stage.)

# 18.09.2012 TH and EK 
# 16.07.2013 EK, BAM sorting changed
# 23.04.2013 MK, polished script and added parameter defining the sort order of paired end data
# 23.04.2013 MK, added standed and add NH-flag parameters

# convert bam to sam, sort bam if the data is paired-end
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
filtersamreads.binary <- file.path(chipster.tools.path, "picard-tools", "FilterSamReads.jar")
samtools.view <- paste(samtools.binary, "view alignment.bam")

if(addnh == "yes" && paired == "yes") {
	#List of IDs having mapping quality >= 4 and having proper pairing (0x0002). Uniq command make sure that both pairs must fulfill filtering criteria
	keep_id.command <- paste(samtools.binary, "view -q4 -F4 -f2 alignment.bam | cut -f1 | sort | uniq -d > alignment_id_keep.txt")
	system(keep_id.command)
	filtersamreads.command <- paste("java -Xmx4096m -jar", filtersamreads.binary, "I=alignment.bam FILTER=includeReadList RLF=alignment_id_keep.txt O=alignment_keep.bam VALIDATION_STRINGENCY=LENIENT")
	system(filtersamreads.command) 
	samtools.view <- paste(samtools.binary, "view alignment_keep.bam | perl -p -e 's/\n$/\tNH:i:1\n/g'")
} else if(addnh == "yes" && paired == "no") {
	samtools.view <- paste(samtools.binary, "view -q4 -F4 alignment.bam | perl -p -e 's/\n$/\tNH:i:1\n/g'")
}

# gtf
gtf <- file.path(chipster.tools.path, "genomes", "gtf", organism)

# exoncount
dexseq.binary <- file.path(chipster.tools.path, "dexseq-exoncounts", "dexseq_count.py")
if(paired == "yes") {
	dexseq.binary <- paste(dexseq.binary, ifelse(sort == "pos", paste("-p yes -r pos"), paste("-p yes -r name")), sep =" ")
}
dexseq.command <- paste("python", dexseq.binary, "-s", stranded, gtf, "- exon-counts-out.tsv")

# run
command <- paste(samtools.view, "|", dexseq.command)
system(command)

# separate result file
system("head -n -4 exon-counts-out.tsv > exon-counts.tsv")
system("tail -n 4 exon-counts-out.tsv > exon-counts-info.tsv")

# bring in file to R environment for formating
file <- c("exon-counts.tsv")
dat <- read.table(file, header=F, sep="\t")
names(dat) <- c("id", "count")

# write result table to output
write.table(dat, file="exon-counts.tsv", col.names=T, quote=F, sep="\t", row.names=F)

file <- c("exon-counts-info.tsv")
dat <- read.table(file, header=F, sep="\t")
names(dat) <- c("id", "count")

# write result table to output
write.table(dat, file="exon-counts-info.tsv", col.names=T, quote=F, sep="\t", row.names=F)

# EOF

