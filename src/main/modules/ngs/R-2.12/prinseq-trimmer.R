# TOOL prinseq-trimmer.R: "Trim reads for several criteria" (Trims reads based on given criteria. This tool is based on the PRINSEQ package.)
# INPUT fastqfile: "Input reads set" TYPE GENERIC
# INPUT OPTIONAL matepair_fastqfile: "Input reads mate pair file" TYPE GENERIC
# OUTPUT OPTIONAL trimmed.fastq
# OUTPUT OPTIONAL trimmed_matepair.fastq
# OUTPUT OPTIONAL trimmed.fasta
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
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "no", y: "yes"] DEFAULT y (Write a log file)

# KM 17.1.2012
# EK 7.5.2013 Reorganized parameters
# AMS 06.11.2013 Added support for paired-end reads

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

# Check if two files were given as input, and if so run the python script
# that interlaces the mate pairs into a single file
input_files <- dir()
is_paired_end <- (length(grep("matepair_fastqfile", input_files))>0)
if (is_paired_end) {
	#check if the fastqfomat is hiseq or old
	first_row_1 <- read.table(file="fastqfile", nrow=1, header=FALSE, sep=" ", check.names=FALSE, comment.char="")
	first_row_2 <- read.table(file="matepair_fastqfile", nrow=1, header=FALSE, sep=" ", check.names=FALSE, comment.char="")
	if ( ncol(first_row_1) == 2 ) { fq.hiseq <- "yes" } else { fq.hiseq <- "no" }
	
	
	# check that the input files are truly matepairs by comparing
	# sequence ID, discarding the last character
	name_length <- nchar(as.character(first_row_1[1,1]))
	id_1 <- substr(as.character(first_row_1[1,1]), start=1, stop=name_length-1)
	id_2 <- substr(as.character(first_row_2[1,1]), start=1, stop=name_length-1)
	if (id_1 != id_2) {
		stop("CHIPSTER-NOTE: It appears that the two input files are not matepairs. Please check that the correct input files were selected.")
	}
	
	# figure out which file is the first and second matepair, and issue
	# the python script call accordingly
	if ( fq.hiseq == "no"){		
		mate_number <- substr(as.character(first_row_1[1,1]), start=name_length, stop=name_length)
	} else {
		mate_number <- substr(as.character(first_row_1[1,2]), start=1, stop=1)
	}	
	if (mate_number == "1") {
		binary_python_scripts <- file.path(chipster.module.path, "shell", "match-mate-pairs", "interleave-fastq.py")
		system_command <- paste("python", binary_python_scripts, "fastqfile", "matepair_fastqfile", "interleaved_fastqfile")
		system(system_command)	
		system("echo Executed interleave python script with: > filter.log")
		echo.command <- paste("echo '", system_command, "'>> filter.log")
		system(echo.command)
	} else {
		binary_python_scripts <- file.path(chipster.module.path, "shell", "match-mate-pairs", "interleave-fastq.py")
		system_command <- paste("python", binary_python_scripts, "matepair_fastqfile", "fastqfile", "interleaved_fastqfile")
		system(system_command)	
		system("echo Executed interleave python script with: > filter.log")
		echo.command <- paste("echo '", system_command, "'>> filter.log")
		system(echo.command)
	}	
	
	# remove input files to clear up disk space
	system("rm -f fastqfile")
	system("rm -f matepair_fastqfile")
	system("mv interleaved_fastqfile fastqfile")
	
}

# binary
binary.prinseq <- c(file.path(chipster.tools.path, "prinseq", "prinseq-lite.pl"))

trim.params <- paste(" ")
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

if (input.mode == "fq") {
 trim.command <- paste(binary.prinseq, trim.params, "-fastq fastqfile -out_good trimmed")
}

if (input.mode == "fa") {
	trim.command <- paste(binary.prinseq, trim.params, "-fasta fastqfile -out_good trimmed")
}


if (log.file == "y") {
	system("echo Running PRINSEQ filtering with command: > trim.log")
	echo.command <- paste("echo '", trim.command, "'>> trim.log")
	system(echo.command)
	trim.command <- paste(trim.command, "-verbose 2>> trim.log")
}


system(trim.command)
#Make sure something is in the output
if (input.mode == "fq") {
	system("if [ ! -e  trimmed.fastq ] ; then echo 'Trimming produced an empty trimmed.fastq sequence set' >> trim.log ; echo '' > trimmed.fastq ; fi")
}

if (input.mode == "fa") {
	system("if [ ! -e  trimmed.fasta ] ; then echo 'Trimming produced an empty trimmed.fasta sequence set' >> trim.log ; echo '' > trimmed.fasta ;fi")
}

# If filtering on paired-end data, perform matching of
# mate pairs using python script and then de-interlace them to two files
if (is_paired_end) {
	if ( fq.hiseq == "no"){	
		binary_python_scripts <- file.path(chipster.module.path, "shell", "match-mate-pairs", "match-pairs.py")
	} else{
		binary_python_scripts <- file.path(chipster.module.path, "shell", "match-mate-pairs", "match-pairs-hiseq.py")
	}
	
	system_command <- paste("python", binary_python_scripts, "trimmed.fastq", "matched_fastqfile >> filter.log")
	system(system_command)
	
	system("echo Executed match_pair python script with: >> filter.log")
	echo.command <- paste("echo '", system_command, "'>> filter.log")
	system(echo.command)
	
	# remove input files to clear up disk space
	system("rm -f accepted.fastq")
	
	binary_python_scripts <- file.path(chipster.module.path, "shell", "match-mate-pairs", "deinterleave-fastq.py")
	system_command <- paste("python", binary_python_scripts, "matched_fastqfile", "trimmed.fastq", "trimmed_matepair.fastq")
	system(system_command)	
	
	system("echo Executed deinterleave python script with: >> filter.log")
	echo.command <- paste("echo '", system_command, "'>> filter.log")
	system(echo.command)
	
	# remove input files to clear up disk space
	system("rm -f matched_fastqfile")
}

#stop(paste('CHIPSTER-NOTE: ', filter.command))

