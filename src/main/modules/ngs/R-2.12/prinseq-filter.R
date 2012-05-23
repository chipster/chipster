# TOOL prinseq-filter.R: "Filter reads for several criteria" (Filters reads based on several criteria. Different criterias are combined with AND operator. This tool is based on the PRINSEQ package.)
# INPUT fastqfile: "Input reads file" TYPE GENERIC
# INPUT OPTIONAL matepair_fastqfile: "Input reads mate pair file" TYPE GENERIC
# OUTPUT OPTIONAL accepted.fastq
# OUTPUT OPTIONAL accepted.fasta
# OUTPUT OPTIONAL accepted_matepair.fastq
# OUTPUT OPTIONAL rejected.fastq
# OUTPUT OPTIONAL rejected.fasta
# OUTPUT OPTIONAL filter.log
# PARAMETER OPTIONAL output.mode: "Results to write out" TYPE [ filt: "accepted sequences only", both: "accepted and rejected sequences into separate files"] DEFAULT filt (With this section you can define if the sequences that get filtered out are collected to a separate file) 
# PARAMETER OPTIONAL input.mode: "Input file format" TYPE [ fq: "FASTQ", fa: "FASTA"] DEFAULT fq (Define the file format of the reads file)
# PARAMETER OPTIONAL noniupac: "Remove reads with non-standard characters" TYPE [yes, no] DEFAULT no (Filter out reads with characters other than A, C, G, T or N.)
# PARAMETER OPTIONAL max.len: "Maximum length" TYPE INTEGER (Select only reads that are shorter than the given value.)
# PARAMETER OPTIONAL min.len: "Minimum length" TYPE INTEGER (Select only reads that are longer than the given value.)
# PARAMETER OPTIONAL max.gc: "Maximum GC content" TYPE INTEGER (Select only reads that have GC content smaller than the given value.)
# PARAMETER OPTIONAL min.gc: "Minimum GC content" TYPE INTEGER (Select only reads that have GC content larger than the given value.)
# PARAMETER OPTIONAL min.qual.score: "Minimum quality score" TYPE INTEGER (Filter reads with GC content below than the given value.)
# PARAMETER OPTIONAL max.qual.score: "Maximum quality score" TYPE INTEGER (Filter reads with GC content above then the given value.)
# PARAMETER OPTIONAL min.qual.mean: "Minimum mean quality" TYPE INTEGER ( Filter reads with quality score mean below the given value.)
# PARAMETER OPTIONAL max.qual.mean: "Maximum mean quality" TYPE INTEGER ( Filter reads with quality score mean above the given value.)
# PARAMETER OPTIONAL ns.max.p: "Maximum percentage of Ns" TYPE INTEGER (Filter reads for which the percentage of Ns is higher than the given value.)
# PARAMETER OPTIONAL ns.max.n: "Maximun count of Ns" TYPE INTEGER (Filter reads which have more Ns than the given value.)
# PARAMETER OPTIONAL seq.num: "Maximum number of reads" TYPE INTEGER (Only keep the given number number of reads that pass all other filters.)
# PARAMETER OPTIONAL derep: "Type of duplicates to filter" TYPE [0: "none", 1: "exact duplicate", 2: "5-prime duplicate", 3: "3-prime duplicate", 4: "reverse complement exact duplicate", 5:"reverse complement 5-prime/3-prime duplicate"] DEFAULT 0 (Type of duplicates to filter.)
# PARAMETER OPTIONAL derep.min: "Number of allowed duplicates" TYPE INTEGER (This option specifies the number of allowed duplicates. For example, to remove reads that occur more than 5 times, you would specify value 6.)
# PARAMETER OPTIONAL lc.dust: "DUST filter threshold" TYPE INTEGER (Use DUST method with the given maximum allowed score, between 0 and 100.)
# PARAMETER OPTIONAL lc.entropy: "ENTROPY filter threshold" TYPE INTEGER (Use ENTROPY method with the given minimum entropy value, between 0 and 100.)
# PARAMETER OPTIONAL phred64: "Base quality encoding" TYPE [ n: "Sanger", y: "Illumina v1.3-1.5"] DEFAULT n (Select \"Sanger" for Illumina v1.8+, Sanger, Roche/454, Ion Torrent and PacBio data.)
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "no", y: "yes"] DEFAULT y (Write a log file)

# Filter fastq and fasta files based on a number of criteria
# KM, EK, 16-04-2012
# MG, 18-04-2012, added matepair functionality
# KM, 22-05-2012, fixed fastq checking

# Check out if the files are compressed and if so unzip it
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
	
	# check that the input files are truly matepairs by comparing
	# sequence ID, discarding the last character
	first_row_1 <- read.table(file="fastqfile", nrow=1, header=FALSE, sep="\t", check.names=FALSE, comment.char="")
	first_row_2 <- read.table(file="matepair_fastqfile", nrow=1, header=FALSE, sep="\t", check.names=FALSE, comment.char="")
	name_length <- nchar(as.character(first_row_1[1,1]))
	id_1 <- substr(as.character(first_row_1[1,1]), start=1, stop=name_length-1)
	id_2 <- substr(as.character(first_row_2[1,1]), start=1, stop=name_length-1)
	if (id_1 != id_2) {
		stop("CHIPSTER-NOTE: It appears that the two input files are not matepairs. Please check that the correct input files were selected.")
	}
	
	# figure out which file is the first and second matepair, and issue
	# the python script call accordingly
	mate_number <- substr(as.character(first_row_1[1,1]), start=name_length, stop=name_length)
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
binary.prinseq <- c(file.path(chipster.tools.path, "prinseq", "prinseq-lite.pl" ))

filter.params <- paste(" ")
if (!is.na(max.len)) {
	filter.params <- paste(filter.params, "-max_len",  max.len)
}

if (!is.na(min.len)) {
	filter.params <- paste(filter.params, "-min_len",  min.len)
}

if (!is.na(max.gc)) {
	filter.params <- paste(filter.params, "-max_gc",  max.gc)
}

if (!is.na(min.gc)) {
	filter.params <- paste(filter.params, "-min_gc",  min.gc)
}

if (!is.na(min.qual.score)) {
	filter.params <- paste(filter.params, "-min_qual_score",  min.qual.score)
}

if (!is.na(max.qual.score)) {
	filter.params <- paste(filter.params, "-max_qual_score",  max.qual.score)
}

if (!is.na(min.qual.mean)) {
	filter.params <- paste(filter.params, "-min_qual_mean",  min.qual.mean)
}

if (!is.na(max.qual.mean)) {
	filter.params <- paste(filter.params, "-max_qual_mean",  max.qual.mean)
}

if (!is.na(ns.max.p)) {
	filter.params <- paste(filter.params, "-ns_max_p",  ns.max.p)
}

if (!is.na(ns.max.n)) {
	filter.params <- paste(filter.params, "-ns_max_n",  ns.max.n)
}

if (!is.na(seq.num)) {
	filter.params <- paste(filter.params, "-seq_num",  seq.num)
}

if (!is.na(lc.dust)) {
	filter.params <- paste(filter.params, "-lc_method dust -lc_threshold", lc.dust )
}

if (!is.na(lc.entropy)) {
	filter.params <- paste(filter.params, "-lc_method entropy -lc_threshold", lc.entropy )
}

if (derep > 0) {	
	filter.params <- paste(filter.params, "-derep", derep )
    if (derep == "1") {
	   filter.params <- paste(filter.params, "-derep_min", derep.min )
    }
    if (derep == "4") {
	   filter.params <- paste(filter.params, "-derep_min", derep.min )
    }
}

if (noniupac == "yes") {
	filter.params <- paste(filter.params, "-noniupac" )
}

if (phred64 == "y") {
	filter.params <- paste(filter.params, "-phred64")
}

if (output.mode == "both") {
	filter.params <- paste(filter.params, "-out_bad rejected")
}

if (input.mode == "fq") {
	filter.command <- paste(binary.prinseq, filter.params, "-fastq fastqfile -out_good accepted")
}

if (input.mode == "fa") {
	filter.command <- paste(binary.prinseq, filter.params, "-fasta fastqfile -out_good accepted")
}

if (log.file == "y") {
	system("echo Running PRINSEQ filtering with command: >> filter.log")
	echo.command <- paste("echo '", filter.command, "'>> filter.log")
	system(echo.command)
	filter.command <- paste(filter.command, "-verbose 2>> filter.log")
}

system(filter.command)

# Make sure something is in the output
if (input.mode == "fq") {
	system("if [ ! -e  accepted.fastq ] ; then echo 'Filtering produced an empty accepted.fastq sequence set' > accepted.fastq ; fi")
}

if (input.mode == "fa") {
	system("if [ ! -e  accepted.fasta ] ; then echo 'Filtering produced an empty accepted.fasta sequence set' > accepted.fasta ; fi")
}

if (output.mode == "both") {
	if (input.mode == "fq") {
		system("if [ ! -e  rejected.fastq ] ; then echo 'Filtering produced an empty rejected.fastq sequence set' >> filter.log ; echo '' > rejected.fastq ; fi")
	}
	
	if (input.mode == "fa") {
		system("if [ ! -e  rejected.fasta ] ; then echo 'Filtering produced an empty rejected.fasta sequence set' >> filter.log  ; echo '' > rejected.fasta ; fi")
	}
}

# remove input files to clear up disk space
system("rm -f fastqfile")

# If filtering on paired-end data, perform matching of
# mate pairs using python script and then de-interlace them to two files
if (is_paired_end) {
	binary_python_scripts <- file.path(chipster.module.path, "shell", "match-mate-pairs", "match-pairs.py")
	system_command <- paste("python", binary_python_scripts, "accepted.fastq", "matched_fastqfile")
	system(system_command)
	
	system("echo Executed match_pair python script with: >> filter.log")
	echo.command <- paste("echo '", system_command, "'>> filter.log")
	system(echo.command)

	# remove input files to clear up disk space
	system("rm -f accepted.fastq")

	binary_python_scripts <- file.path(chipster.module.path, "shell", "match-mate-pairs", "deinterleave-fastq.py")
	system_command <- paste("python", binary_python_scripts, "matched_fastqfile", "accepted.fastq", "accepted_matepair.fastq")
	system(system_command)	

	system("echo Executed deinterleave python script with: >> filter.log")
	echo.command <- paste("echo '", system_command, "'>> filter.log")
	system(echo.command)

	# remove input files to clear up disk space
	system("rm -f matched_fastqfile")
}

# stop
