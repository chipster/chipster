# TOOL tophat2-single-end-with-index-building.R: "TopHat2 for single end reads and own genome" (Aligns Illumina RNA-seq reads to a genome in order to identify exon-exon splice junctions. The alignment process consists of several steps. If annotation is available as a GTF file, TopHat will extract the transcript sequences and use Bowtie to align reads to this virtual transcriptome first. Only the reads that do not fully map to the transcriptome will then be mapped on the genome. The reads that still remain unmapped are split into shorter segments, which are then aligned to the genome. Alignment results are given in a BAM file, which is automatically indexed and hence ready to be viewed in Chipster genome browser.)
# INPUT reads1.fq: "Reads" TYPE GENERIC
# INPUT genome.txt: "Genome to align against" TYPE GENERIC
# INPUT OPTIONAL genes.gtf: "Optional GTF file" TYPE GENERIC
# OUTPUT OPTIONAL tophat.bam
# OUTPUT OPTIONAL tophat.bam.bai
# OUTPUT OPTIONAL junctions.bed
# OUTPUT OPTIONAL tophat-summary.txt
# OUTPUT OPTIONAL tophat2.log
# PARAMETER OPTIONAL no.novel.juncs: "When GTF file is used, ignore novel junctions" TYPE [yes, no] DEFAULT no (Only look for reads across junctions indicated in the supplied GTF file.)
# PARAMETER OPTIONAL quality.format: "Base quality encoding used" TYPE [sanger: "Sanger - Phred+33", phred64: "Phred+64"] DEFAULT sanger (Quality encoding used in the fastq file.)
# PARAMETER OPTIONAL max.multihits: "How many hits is a read allowed to have" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 20 (Instructs TopHat to allow up to this many alignments to the reference for a given read.)
# PARAMETER OPTIONAL mismatches: "Number of mismatches allowed in final alignment" TYPE INTEGER FROM 0 TO 5 DEFAULT 2 (Final read alignments having more than this many mismatches are discarded.)
# PARAMETER OPTIONAL min.anchor.length: "Minimum anchor length" TYPE INTEGER FROM 3 TO 1000 DEFAULT 8 (TopHat2 will report junctions spanned by reads with at least this many bases on each side of the junction. Note that individual spliced alignments may span a junction with fewer than this many bases on one side. However, every junction involved in spliced alignments is supported by at least one read with this many bases on each side.)
# PARAMETER OPTIONAL splice.mismatches: "Maximum number of mismatches allowed in the anchor" TYPE INTEGER FROM 0 TO 2 DEFAULT 0 (The maximum number of mismatches that may appear in the anchor region of a spliced alignment.)
# PARAMETER OPTIONAL min.intron.length: "Minimum intron length" TYPE INTEGER FROM 10 TO 1000 DEFAULT 70 (TopHat2 will ignore donor-acceptor pairs closer than this many bases apart.)
# PARAMETER OPTIONAL max.intron.length: "Maximum intron length" TYPE INTEGER FROM 1000 TO 1000000 DEFAULT 500000 (TopHat2 will ignore donor-acceptor pairs farther than this many bases apart, except when such a pair is supported by a split segment alignment of a long read.)
# PARAMETER OPTIONAL library.type: "Library type" TYPE [fr-unstranded: fr-unstranded, fr-firststrand: fr-firststrand, fr-secondstrand: fr-secondstrand] DEFAULT fr-unstranded (Which library type to use.)

# EK 17.4.2012 added -G and -g options
# MG 24.4.2012 added ability to use gtf files from Chipster server
# AMS 19.6.2012 added unzipping
# AMS 27.6.2012 added parameter mate.std.dev, allow negative values for mate.inner.distance
# AMS 4.10.2012 added BED sorting
# KM 10.7. 2012 added RN5
# AMS 11.11.2013 added thread support
# AMS 3.1.2014 added transcriptome index for human
# EK 3.1.2014 added alignment summary to output, added quality and mismatch parameter
# AMS 22.5.2014 modified to use own genome
# ML 15.01.2015 Added the library-type parameter
# AMS 29.01.2015 Removed optional outputs deletions.bed and insertions.bed

# OUTPUT OPTIONAL tophat2.log

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads1.fq")
unzipIfGZipFile("genome.txt")
unzipIfGZipFile("genes.gtf")

options(scipen = 10)
# max.intron.length <- formatC(max.intron.length, "f", digits = 0)

# setting up TopHat
tophat.binary <- c(file.path(chipster.tools.path, "tophat2", "tophat2"))
bowtie.binary <- c(file.path(chipster.tools.path, "bowtie2", "bowtie2"))
bowtie2.index.binary <- file.path(chipster.module.path, "shell", "check_bowtie2_index.sh")
path.bowtie <- c(file.path(chipster.tools.path, "bowtie2"))
path.samtools <- c(file.path(chipster.tools.path, "samtools"))
set.path <-paste(sep="", "PATH=", path.bowtie, ":", path.samtools, ":$PATH")

# Do indexing
print("Indexing the genome...")
system("echo Indexing the genome... > bowtie2.log")
check.command <- paste ( bowtie2.index.binary, "genome.txt| tail -1 ")
genome.dir <- system(check.command, intern = TRUE)
bowtie2.genome <- file.path( genome.dir , "genome.txt")


# command start
command.start <- paste("bash -c '", set.path, tophat.binary)

# parameters
#command.parameters <- paste("--bowtie1 -r", mate.inner.distance, "--mate-std-dev", mate.std.dev, "-a", min.anchor.length, "-m", splice.mismatches, "-i", min.intron.length, "-I", max.intron.length, "-g", max.multihits, "--library-type fr-unstranded")
# command.parameters <- paste("-p", chipster.threads.max, "--read-mismatches", mismatches, "-a", min.anchor.length, "-m", splice.mismatches, "-i", min.intron.length, "-I", max.intron.length, "-g", max.multihits, "--library-type fr-unstranded")
command.parameters <- paste("-p", chipster.threads.max, "--read-mismatches", mismatches, "-a", min.anchor.length, "-m", splice.mismatches, "-i", min.intron.length, "-I", max.intron.length, "-g", max.multihits)

if ( quality.format == "phred64") {
	command.parameters <- paste(command.parameters, "--phred64-quals")
}

# optional GTF command, if a GTF file has been provided by user
if (file.exists("genes.gtf")){
	if (no.novel.juncs == "yes") {
		command.parameters <- paste(command.parameters, "-G", "genes.gtf", "--no-novel-juncs")
	} else {
		command.parameters <- paste(command.parameters, "-G", "genes.gtf")
	}
}


if (library.type == "fr-unstranded") {
	command.parameters <- paste(command.parameters, "--library-type fr-unstranded")
}else if (library.type == "fr-firststrand") {
	command.parameters <- paste(command.parameters, "--library-type fr-firststrand")	
}else if (library.type == "fr-secondstrand") {	
	command.parameters <- paste(command.parameters, "--library-type fr-secondstrand")
}

# command ending
#command.end <- paste(path.bowtie.index, "reads1.fq reads2.fq >> tophat2.log '")
command.end <- paste(bowtie2.genome, "reads1.fq 2>>tophat.log'")

# run tophat
command <- paste(command.start, command.parameters, command.end)

echo.command <- paste("echo '",command ,"' 2>> tophat.log " )
system(echo.command)
system("echo >> tophat.log")
#stop(paste('CHIPSTER-NOTE: ', command))
system(command)


# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# sort bam (removed because TopHat itself does the sorting)
# system(paste(samtools.binary, "sort tophat_out/accepted_hits.bam tophat"))
system("mv tophat_out/accepted_hits.bam tophat.bam") 

# index bam
system(paste(samtools.binary, "index tophat.bam"))

system("mv tophat_out/junctions.bed junctions.u.bed")
system("mv tophat_out/insertions.bed insertions.u.bed")
system("mv tophat_out/deletions.bed deletions.u.bed")
system("mv tophat_out/align_summary.txt tophat-summary.txt")

# sorting BEDs
source(file.path(chipster.common.path, "bed-utils.R"))


if (file.exists("junctions.u.bed")){
	size <- file.info("junctions.u.bed")$size
	if (size > 100){	
		bed <- read.table(file="junctions.u.bed", skip=1, sep="\t")
		colnames(bed)[1:2] <- c("chr", "start")
		sorted.bed <- sort.bed(bed)
		write.table(sorted.bed, file="junctions.bed", sep="\t", row.names=F, col.names=F, quote=F)
	}	
}

if (file.exists("insertions.u.bed")){
	size <- file.info("insertions.u.bed")$size
	if (size > 100){
		bed <- read.table(file="insertions.u.bed", skip=1, sep="\t")
		colnames(bed)[1:2] <- c("chr", "start")
		sorted.bed <- sort.bed(bed)
		write.table(sorted.bed, file="insertions.bed", sep="\t", row.names=F, col.names=F, quote=F)
	}
}

if (file.exists("deletions.u.bed")){
	size <- file.info("deletions.u.bed")$size
	if (size > 100){
		bed <- read.table(file="deletions.u.bed", skip=1, sep="\t")
		colnames(bed)[1:2] <- c("chr", "start")
		sorted.bed <- sort.bed(bed)
		write.table(sorted.bed, file="deletions.bed", sep="\t", row.names=F, col.names=F, quote=F)
	}
}

if (!(file.exists("tophat-summary.txt"))){
	#system("mv tophat_out/logs/tophat.log tophat2.log")
	system("mv tophat.log tophat2.log")
}