# TOOL tophat2.R: "TopHat2 for paired end reads" (Aligns paired end Illumina RNA-Seq reads to a genome in order to identify exon-exon splice junctions. The alignment process consists of three steps. If annotation is available as a GTF file, TopHat will extract the transcript sequences and use Bowtie to align reads to this virtual transcriptome first. Only the reads that do not fully map to the transcriptome will then be mapped on the genome. The reads that remain still unmapped are split into shorter segments, which are then aligned to the genome. Segment mappings are used to find potential splice sites. Sequences flanking a splice site are concatenated, and unmapped segments are mapped to them. Segment alignments are then stitched together to form whole read alignments. Alignment results are given in a BAM file, which is automatically indexed and hence ready to be viewed in Chipster genome browser.)
# INPUT reads1.fq: "Reads to align" TYPE GENERIC
# INPUT reads2.fq: "Reads to align" TYPE GENERIC
# INPUT OPTIONAL genes.gtf: "Optional GTF file" TYPE GENERIC
# OUTPUT tophat.bam
# OUTPUT tophat.bam.bai
# OUTPUT OPTIONAL junctions.bed
# OUTPUT OPTIONAL insertions.bed
# OUTPUT OPTIONAL deletions.bed
# OUTPUT OPTIONAL tophat2.log
# PARAMETER genome: "Genome" TYPE [hg19: "Human genome (hg19\)", mm10: "Mouse genome (mm10\)", Rattus_norvegicus.Rnor_5.0.70.dna.toplevel: "Rat genome (rn5\)", rn4: "Rat genome (rn4\)", Halorubrum_lacusprofundi_ATCC_49239: "Halorubrum lacusprofundi ATCC 49239 genome", Canis_familiaris.CanFam3.1.71.dna.toplevel: "Dog genome (Ensembl canFam3\)", Sus_scrofa.Sscrofa10.2.69.dna.toplevel: "Pig (sus_scrofa10.2.69\)",  Gasterosteus_aculeatus.BROADS1.71.dna.toplevel: "Gasterosteus aculeatus genome (BROADS1.71\)", athaliana.TAIR10: "A. thaliana genome (TAIR10\)", Arabidopsis_lyrata.v.1.0.16: "A. lyrata genome (1.0.16\)", ovis_aries_texel: "Sheep genome (oar3.1\)"] DEFAULT hg19 (Genome or transcriptome that you would like to align your reads against.)
# PARAMETER OPTIONAL use.gtf: "Use annotation GTF" TYPE [yes, no] DEFAULT yes (If this option is provided, TopHat will extract the transcript sequences and use Bowtie to align reads to this virtual transcriptome first. Only the reads that do not fully map to the transcriptome will then be mapped on the genome. The reads that did map on the transcriptome will be converted to genomic mappings (spliced as needed\) and merged with the novel mappings and junctions in the final TopHat output. If no GTF file is provided by user, internal annotation file will be used if available. Internal annotation is provided for human, mouse and rat.)
# PARAMETER OPTIONAL no.novel.juncs: "When GTF file is used, ignore novel junctions" TYPE [yes, no] DEFAULT yes (If annotation GTF is used, TopHat will extract the transcript sequences and use Bowtie to align reads to this virtual transcriptome first. Only the reads that do not fully map to the transcriptome will then be mapped on the genome. The reads that did map on the transcriptome will be converted to genomic mappings (spliced as needed\) and merged with the novel mappings and junctions in the final TopHat output. If no GTF file is provided by the user, internal annotation file will be used if available. Internal annotation is provided for human, mouse and rat.)
# PARAMETER OPTIONAL mate.inner.distance: "Expected inner distance between mate pairs" TYPE INTEGER DEFAULT 200 (Expected mean inner distance between mate pairs. For example, if your fragment size is 300 bp and read length is 50 bp, the inner distance is 200.)
# PARAMETER OPTIONAL mate.std.dev: "Standard deviation for the inner distances between mate pairs" TYPE INTEGER DEFAULT 20 (The standard deviation for the distribution on inner distances between mate pairs. The default is 20bp.)
# PARAMETER OPTIONAL min.anchor.length: "Minimum anchor length" TYPE INTEGER FROM 3 TO 1000 DEFAULT 8 (TopHat2 will report junctions spanned by reads with at least this many bases on each side of the junction. Note that individual spliced alignments may span a junction with fewer than this many bases on one side. However, every junction involved in spliced alignments is supported by at least one read with this many bases on each side.)
# PARAMETER OPTIONAL splice.mismatches: "Maximum number of mismatches allowed in the anchor" TYPE INTEGER FROM 0 TO 2 DEFAULT 0 (The maximum number of mismatches that may appear in the anchor region of a spliced alignment.)
# PARAMETER OPTIONAL min.intron.length: "Minimum intron length" TYPE INTEGER FROM 10 TO 1000 DEFAULT 70 (TopHat2 will ignore donor-acceptor pairs closer than this many bases apart.)
# PARAMETER OPTIONAL max.intron.length: "Maximum intron length" TYPE INTEGER FROM 1000 TO 1000000 DEFAULT 500000 (TopHat2 will ignore donor-acceptor pairs farther than this many bases apart, except when such a pair is supported by a split segment alignment of a long read.)
# PARAMETER OPTIONAL max.multihits: "How many hits is a read allowed to have" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 20 (Instructs TopHat to allow up to this many alignments to the reference for a given read, and suppresses all alignments for reads with more than this many alignments.)

# EK 17.4.2012 added -G and -g options
# MG 24.4.2012 added ability to use gtf files from Chipster server
# AMS 19.6.2012 added unzipping
# AMS 27.6.2012 added parameter mate.std.dev, allow negative values for mate.inner.distance
# AMS 4.10.2012 added BED sorting
# KM 10.7. 2012 Added RN5


# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads1.fq")
unzipIfGZipFile("reads2.fq")

options(scipen = 10)
# max.intron.length <- formatC(max.intron.length, "f", digits = 0)

# setting up TopHat
tophat.binary <- c(file.path(chipster.tools.path, "tophat2", "tophat2"))
path.bowtie <- c(file.path(chipster.tools.path, "bowtie2"))
path.samtools <- c(file.path(chipster.tools.path, "samtools"))
set.path <-paste(sep="", "PATH=", path.bowtie, ":", path.samtools, ":$PATH")
path.bowtie.index <- c(file.path(path.bowtie, "indexes", genome))

# command start
command.start <- paste("bash -c '", set.path, tophat.binary)

# parameters
#command.parameters <- paste("--bowtie1 -r", mate.inner.distance, "--mate-std-dev", mate.std.dev, "-a", min.anchor.length, "-m", splice.mismatches, "-i", min.intron.length, "-I", max.intron.length, "-g", max.multihits, "--library-type fr-unstranded")
command.parameters <- paste("-r", mate.inner.distance, "--mate-std-dev", mate.std.dev, "-a", min.anchor.length, "-m", splice.mismatches, "-i", min.intron.length, "-I", max.intron.length, "-g", max.multihits, "--library-type fr-unstranded")

# optional GTF command, if a GTF file has been provided by user
command.gtf <- ""
input_files <- dir()
is_gtf <- (length(grep("genes.gtf", input_files))>0)
if (is_gtf) {
	if (no.novel.juncs == "yes") {
		command.gtf <- paste("-G", "genes.gtf", "--no-novel-juncs")
	} else {
		command.gtf <- paste("-G", "genes.gtf")
	}
}

# optional GTF command, if a GTF file has NOT been provided by user
# BUT is avaliable from Chipster server
genome_available <- FALSE
if (genome == "hg19" ||	genome == "mm9" || genome == "mm10" || genome == "rn4" || genome == "Rattus_norvegicus.Rnor_5.0.70.dna.toplevel" || genome == "Canis_familiaris.CanFam3.1.71.dna.toplevel") genome_available <- TRUE
if (!is_gtf && genome_available) {

	# annotation file setup
	if (genome == "hg19") {
		annotation.file <- "Homo_sapiens.GRCh37.68.chr.gtf"
	}
	if (genome == "mm10") {
		annotation.file <- "Mus_musculus.GRCm38.68.chr.gtf"
	}
	if (genome == "mm9") {
		annotation.file <- "Mus_musculus.NCBIM37.62.chr.gtf"
	}
	if (genome == "rn4") {
		annotation.file <- "Rattus_norvegicus.RGSC3.4.62.chr.gtf"
	}
	if (genome == "Rattus_norvegicus.Rnor_5.0.70.dna.toplevel") {
		annotation.file <- "Rattus_norvegicus.Rnor_5.0.70.gtf"
	}
	if (genome == "Canis_familiaris.CanFam3.1.71.dna.toplevel") {
		annotation.file <- "Canis_familiaris.CanFam3.1.71.gtf"
	}
	annotation.file <- c(file.path(chipster.tools.path, "genomes", "gtf", annotation.file))	
	
	if (no.novel.juncs == "yes") {
		command.gtf <- paste("-G", annotation.file, "--no-novel-juncs")
	} else {
		command.gtf <- paste("-G", annotation.file)
	}
}

# command ending
command.end <- paste(path.bowtie.index, "reads1.fq reads2.fq >> tophat2.log '")

# run tophat
if (use.gtf == "yes"){ 
	command <- paste(command.start, command.parameters, command.gtf, command.end)
}else{
	command <- paste(command.start, command.parameters, command.end)
}
echo.command <- paste("echo '",command ,"' > tophat2.log " )
system(echo.command)
system(command)

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# sort bam
system(paste(samtools.binary, "sort tophat_out/accepted_hits.bam tophat"))

# index bam
system(paste(samtools.binary, "index tophat.bam"))

system("mv tophat_out/junctions.bed junctions.u.bed")
system("mv tophat_out/insertions.bed insertions.u.bed")
system("mv tophat_out/deletions.bed deletions.u.bed")

# sorting BEDs
source(file.path(chipster.common.path, "bed-utils.R"))

size <- 0
size <- file.info("junctions.u.bed")$size

echo.command <- paste("echo ",size ," >> tophat2.log" )
system(echo.command)
#system("sleep 120")
#stop(paste('CHIPSTER-NOTE: ', size))



if (size > 100){	
	bed <- read.table(file="junctions.u.bed", skip=1, sep="\t")
	colnames(bed)[1:2] <- c("chr", "start")
	sorted.bed <- sort.bed(bed)
	write.table(sorted.bed, file="junctions.bed", sep="\t", row.names=F, col.names=F, quote=F)
}

size <- file.info("insertions.u.bed")$size
if (size > 100){
	bed <- read.table(file="insertions.u.bed", skip=1, sep="\t")
	colnames(bed)[1:2] <- c("chr", "start")
	sorted.bed <- sort.bed(bed)
	write.table(sorted.bed, file="insertions.bed", sep="\t", row.names=F, col.names=F, quote=F)
}

size <- file.info("deletions.u.bed")$size
if (size > 100){
	bed <- read.table(file="deletions.u.bed", skip=1, sep="\t")
	colnames(bed)[1:2] <- c("chr", "start")
	sorted.bed <- sort.bed(bed)
	write.table(sorted.bed, file="deletions.bed", sep="\t", row.names=F, col.names=F, quote=F)
}
