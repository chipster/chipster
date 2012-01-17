# TOOL tophat.R: "TopHat for paired end reads" (Aligns paired end Illumina RNA-Seq reads to a genome in order to identify exon-exon splice junctions. TopHat first identifies potential exons by mapping the reads to the genome. Using this initial mapping, it builds a database of possible splice junctions, and then maps the reads against these junctions to confirm them. The alignment results are given in a BAM file which is automatically indexed and hence ready to be viewed in Chipster genome browser.)
# INPUT reads1.fq: "Reads to align" TYPE GENERIC
# INPUT reads2.fq: "Reads to align" TYPE GENERIC
# OUTPUT tophat.bam
# OUTPUT tophat.bam.bai
# OUTPUT junctions.bed
# OUTPUT insertions.bed
# OUTPUT deletions.bed
# PARAMETER genome: "Genome" TYPE [hg19: "Human (hg19\)", mm9: "Mouse (mm9\)", rn4: "Rat (rn4\)", Phytophthora_infestans1_1.12: "Phytophthora infestans genome 1.1.12", Populus_trichocarpa.JGI2.0.12: "Populus trichocarpa genome JGI2.0.12"] DEFAULT mm9 (Genome that you would like to align your reads against.)
# PARAMETER mate.inner.distance: "Expected inner distance between mate pairs" TYPE INTEGER FROM 10 TO 1000 DEFAULT 200 (Expected mean inner distance between mate pairs. For example, if your fragment size is 300 bp and read length is 50 bp, the inner distance is 200.)
# PARAMETER OPTIONAL min.anchor.length: "Minimum anchor length" TYPE INTEGER FROM 3 TO 1000 DEFAULT 8 (TopHat will report junctions spanned by reads with at least this many bases on each side of the junction. Note that individual spliced alignments may span a junction with fewer than this many bases on one side. However, every junction involved in spliced alignments is supported by at least one read with this many bases on each side.)
# PARAMETER OPTIONAL splice.mismatches: "Maximum number of mismatches allowed in the anchor" TYPE INTEGER FROM 0 TO 2 DEFAULT 0 (The maximum number of mismatches that may appear in the anchor region of a spliced alignment.)
# PARAMETER OPTIONAL min.intron.length: "Minimum intron length" TYPE INTEGER FROM 10 TO 1000 DEFAULT 70 (TopHat will ignore donor-acceptor pairs closer than this many bases apart.)
# PARAMETER OPTIONAL max.intron.length: "Maximum intron length" TYPE INTEGER FROM 1000 TO 1000000 DEFAULT 500000 (TopHat will ignore donor-acceptor pairs farther than this many bases apart, except when such a pair is supported by a split segment alignment of a long read.)
# PARAMETER OPTIONAL min.isoform.fraction: "Minimum isoform fraction" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.15 (TopHat filters out junctions supported by too few alignments. Suppose a junction spanning two exons is supported by S reads. Let the average depth of coverage of exon A be D, and assume that it is higher than B. If S divided by D is less than the minimum isoform fraction, the junction is not reported. A value of zero disables the filter.)

options(scipen = 10)
# max.intron.length <- formatC(max.intron.length, "f", digits = 0)

# setting up TopHat
tophat.binary <- c(file.path(chipster.tools.path, "tophat", "tophat"))
path.bowtie <- c(file.path(chipster.tools.path, "bowtie"))
path.samtools <- c(file.path(chipster.tools.path, "samtools"))
set.path <-paste(sep="", "PATH=", path.bowtie, ":", path.samtools, ":$PATH")
path.bowtie.index <- c(file.path(path.bowtie, "indexes", genome))


# command start
command.start <- paste("bash -c '", set.path, tophat.binary)

# parameters
command.parameters <- paste("-r", mate.inner.distance, "-a", min.anchor.length, "-m", splice.mismatches, "-i", min.intron.length, "-I", max.intron.length, "-F", min.isoform.fraction, "--library-type fr-unstranded")


# command ending
command.end <- paste(path.bowtie.index, "reads1.fq reads2.fq'")

# run tophat
command <- paste(command.start, command.parameters, command.end)
system(command)

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# sort bam
system(paste(samtools.binary, "sort tophat_out/accepted_hits.bam tophat"))

# index bam
system(paste(samtools.binary, "index tophat.bam"))

system("mv tophat_out/junctions.bed junctions.bed")
system("mv tophat_out/insertions.bed insertions.bed")
system("mv tophat_out/deletions.bed deletions.bed")

