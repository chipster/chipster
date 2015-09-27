# TOOL prinseq-trimmer.R: "Trim reads for several criteria with PRINSEQ" (Trims reads based on given criteria. This tool is based on the PRINSEQ package.)
# INPUT fastqfile: "Input reads set" TYPE GENERIC
# INPUT OPTIONAL matepair_fastqfile: "Input reads mate pair file" TYPE GENERIC
# OUTPUT OPTIONAL trimmed.fastq.gz
# OUTPUT OPTIONAL trimmed_1.fastq.gz
# OUTPUT OPTIONAL trimmed_2.fastq.gz
# OUTPUT OPTIONAL trimmed_1_singletons.fastq.gz
# OUTPUT OPTIONAL trimmed_2_singletons.fastq.gz
# OUTPUT OPTIONAL trimmed.fasta.gz
# OUTPUT OPTIONAL trimmed_1.fasta.gz
# OUTPUT OPTIONAL trimmed_2.fasta.gz
# OUTPUT OPTIONAL trimmed_1_singletons.fasta.gz
# OUTPUT OPTIONAL trimmed_2_singletons.fasta.gz
# OUTPUT OPTIONAL trim.log
# PARAMETER input.mode: "Input file format" TYPE [ fq: "FASTQ", fa: "FASTA"] DEFAULT fq (Define the file format of the reads file)
# PARAMETER phred64: "Base quality encoding" TYPE [ n: "Sanger", y: "Illumina v1.3-1.5"] DEFAULT n (Select \"Sanger" for Illumina v1.8+, Sanger, Roche/454, Ion Torrent and PacBio data.)
# PARAMETER OPTIONAL trim.to.len: "Trim to length" TYPE INTEGER (Trim reads from 3-prime end to this length.)
# PARAMETER OPTIONAL trim.left: "Trim left" TYPE INTEGER (Trim reads at the 5-prime end by given number of bases.)
# PARAMETER OPTIONAL trim.right: "Trim right" TYPE INTEGER (Trim reads at the 3-prime end by given number of bases.)
# PARAMETER OPTIONAL trim.qual.left: "Trim 5-prime end by quality" TYPE INTEGER (Trim reads from the 5-prime end based on the given quality threshold score.)
# PARAMETER OPTIONAL trim.qual.right: "Trim 3-prime end by quality" TYPE INTEGER (Trim reads from the 3-prime end based on the given quality threshold score.)		
# PARAMETER OPTIONAL trim.qual.type: "Quality score calculation method" TYPE [min: "minimum quality value", mean: "mean quality", max: "maximum quality value", sum: "total sum of quality values"] DEFAULT min (Type of quality score calculation to use.)
# PARAMETER OPTIONAL trim.qual.rule: "Quality score comparison condition" TYPE [ lt: "less than", gr: "greater than", et: "equal to"] DEFAULT lt (Rule to use to compare quality score threshold to calculated value.)
# PARAMETER OPTIONAL trim.qual.window: "Window size for quality calculation" TYPE INTEGER DEFAULT 1 (The sliding window size used to calculate quality score by type. Use 1 to stop at the first base that fails the defined rule.)
# PARAMETER OPTIONAL trim.qual.step: "Step size used to move the quality window" TYPE INTEGER DEFAULT 1 (Step size to move the sliding window. To move the window over all quality scores without missing any, the step size should be less than or equal to the window size.)
# PARAMETER OPTIONAL trim.tail.left: "Trim left A/T tails" TYPE INTEGER (Trim poly-A/T tail with a minimum length of the given value at the 5-prime end.)
# PARAMETER OPTIONAL trim.tail.right: "Trim right A/T tails" TYPE INTEGER (Trim poly-A/T tail with a minimum length of the given value at the 3-prime end.)
# PARAMETER OPTIONAL trim.ns.left: "Trim left poly-N tails" TYPE INTEGER (Trim poly-N tail with a minimum length of the given value at the 5-prime end.)
# PARAMETER OPTIONAL trim.ns.right: "Trim right poly-N tails" TYPE INTEGER (Trim poly-N tail with a minimum length of the given value at the 3-prime end.)		
# PARAMETER OPTIONAL min.len: "Minimum length" TYPE INTEGER (Select only reads that are longer than the given value after trimming.)
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "no", y: "yes"] DEFAULT y (Write a log file)
# PARAMETER OPTIONAL singletons: "Write singletons for paired end reads" TYPE [yes, no] DEFAULT no (Write singletons in separate files for paired end reads.)


# KM 17.1.2012
# EK 7.5.2013 Reorganized parameters
# AMS 06.11.2013 Added support for paired-end reads
# AMS 07.01.2014 Paired-end reads now handled by PRINSEQ
# AMS 17.2.2014, gzip outputs

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("fastqfile")
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("matepair_fastqfile")

# Check whether input files are fastq
if (input.mode == "fq") {
	first_four_rows <- read.table(file="fastqfile", nrow=4, header=FALSE, sep="\t", check.names=FALSE, comment.char="")
	## compare sequence ID with quality score id, but discard first character
	#name_length <- nchar(as.character(first_four_rows[1,1]))
	seq_char <- substr(as.character(first_four_rows[1,1]), start=1, stop=1)
	quality_char <- substr(as.character(first_four_rows[3,1]), start=1, stop=1)
	if (seq_char != "@") {
		stop("CHIPSTER-NOTE: It appears as though the input file(s) are not in fastq format. Please check input files or rerun the tool but with the 'Input file format' parameter set to 'FASTA'.")
	}
	if (quality_char != "+") {
		stop("CHIPSTER-NOTE: It appears as though the input file(s) are not in fastq format. Please check input files or rerun the tool but with the 'Input file format' parameter set to 'FASTA'.")
	}
}

# Check if two files were given as input
input_files <- dir()
is_paired_end <- (length(grep("matepair_fastqfile", input_files))>0)

# binary
binary.prinseq <- c(file.path(chipster.tools.path, "prinseq", "prinseq-lite.pl"))

# Parameters
trim.params <- paste("")
if (!is.na(trim.to.len)) {
	trim.params <- paste(trim.params, "-trim_to_len", trim.to.len )
}

if (!is.na(trim.left)) {
	trim.params <- paste(trim.params, "-trim_left", trim.left )
}

if (!is.na(trim.right)) {
	trim.params <- paste(trim.params, "-trim_right", trim.right )
}

if (!is.na(trim.tail.left)) {
	trim.params <- paste(trim.params, "-trim_tail_left",  trim.tail.left)
}

if (!is.na(trim.tail.right)) {
	trim.params <- paste(trim.params, "-trim_tail_right",  trim.tail.right)
}

if (!is.na(trim.ns.left)) {
	trim.params <- paste(trim.params, "-trim_ns_left",  trim.ns.left)
}

if (!is.na(trim.ns.right)) {
	trim.params <- paste(trim.params, "-trim_ns_right",  trim.ns.right)
}

if (!is.na(trim.qual.right)) {
	trim.params <- paste(trim.params, "-trim_qual_right",  trim.qual.right, "-trim_qual_type", trim.qual.type, "-trim_qual_rule", trim.qual.rule, "-trim_qual_window", trim.qual.window, "-trim_qual_step",  trim.qual.step)
}

if (!is.na(trim.qual.left)) {
	trim.params <- paste(trim.params, "-trim_qual_left",  trim.qual.left, "-trim_qual_type", trim.qual.type, "-trim_qual_rule", trim.qual.rule, "-trim_qual_window", trim.qual.window, "-trim_qual_step",  trim.qual.step)
}

if ( phred64 == "y") {
	trim.params <- paste(trim.params, "-phred64")
}

if (!is.na(min.len)) {
	trim.params <- paste(trim.params, "-min_len",  min.len)
}

if (input.mode == "fq") {
	trim.params <- paste(trim.params, "-no_qual_header -fastq fastqfile")
	if (is_paired_end) {
		trim.params <- paste(trim.params, "-fastq2 matepair_fastqfile")
	} 
}

if (input.mode == "fa") {
	trim.params <- paste(trim.params, "-fasta fastqfile")
	if (is_paired_end) {
		trim.params <- paste(trim.params, "-fasta2 matepair_fastqfile")
	}
}
trim.command <- paste(binary.prinseq, trim.params, "-out_good trimmed")

if (log.file == "y") {
	system("echo Running PRINSEQ filtering with command: > trim.log")
	echo.command <- paste("echo '", trim.command, "'>> trim.log")
	system(echo.command)
	trim.command <- paste(trim.command, "-verbose 2>> trim.log")
}


system(trim.command)

# There is no option in PRINSEQ to not write the singletons files, so if they are not required, we delete them.
if (singletons == "no"){
	system("rm -f *_singletons.*")
}

# Compress output files
system("gzip *.fastq")
system("gzip *.fasta")


# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

base1 <- strip_name(inputnames$fastqfile)
if (is_paired_end) {
	base2 <- strip_name(inputnames$matepair_fastqfile)
}else{
	base2 <- ""
}

# Make a matrix of output names
outputnames <- matrix(NA, nrow=10, ncol=2)

# SE fastq
outputnames[1,] <- c("trimmed.fastq.gz", paste(base1, ".fastq.gz", sep =""))
# SE fasta
outputnames[2,] <- c("trimmed.fasta.gz", paste(base1, ".fasta.gz", sep =""))
# PE fastq
outputnames[3,] <- c("trimmed_1.fastq.gz", paste(base1, ".fastq.gz", sep =""))
outputnames[4,] <- c("trimmed_1_singletons.fastq.gz", paste(base1, "_singletons.fastq.gz", sep =""))
outputnames[5,] <- c("trimmed_2.fastq.gz", paste(base2, ".fastq.gz", sep =""))
outputnames[6,] <- c("trimmed_2_singletons.fastq.gz", paste(base2, "_singletons.fastq.gz", sep =""))
# PE fasta
outputnames[7,] <- c("trimmed_1.fasta.gz", paste(base1, ".fasta.gz", sep =""))
outputnames[8,] <- c("trimmed_1_singletons.fasta.gz", paste(base1, "_singletons.fasta.gz", sep =""))
outputnames[9,] <- c("trimmed_2.fasta.gz", paste(base2, ".fasta.gz", sep =""))
outputnames[10,] <- c("trimmed_2_singletons.fasta.gz", paste(base2, "_singletons.fasta.gz", sep =""))

# Write output definitions file
write_output_definitions(outputnames)
