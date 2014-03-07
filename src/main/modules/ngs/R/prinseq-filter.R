# TOOL prinseq-filter.R: "Filter reads for several criteria with PRINSEQ" (Filters reads based on several criteria. Different criterias are combined with AND operator. This tool is based on the PRINSEQ package.)
# INPUT fastqfile: "Input reads file" TYPE GENERIC
# INPUT OPTIONAL matepair_fastqfile: "Input reads mate pair file" TYPE GENERIC
# OUTPUT OPTIONAL accepted.fastq.gz
# OUTPUT OPTIONAL accepted_1.fastq.gz
# OUTPUT OPTIONAL accepted_1_singletons.fastq.gz
# OUTPUT OPTIONAL accepted_2.fastq.gz
# OUTPUT OPTIONAL accepted_2_singletons.fastq.gz
# OUTPUT OPTIONAL accepted.fasta.gz
# OUTPUT OPTIONAL accepted_1.fasta.gz
# OUTPUT OPTIONAL accepted_1_singletons.fasta.gz
# OUTPUT OPTIONAL accepted_2.fasta.gz
# OUTPUT OPTIONAL accepted_2_singletons.fasta.gz
# OUTPUT OPTIONAL rejected.fastq.gz
# OUTPUT OPTIONAL rejected_1.fastq.gz
# OUTPUT OPTIONAL rejected_2.fastq.gz
# OUTPUT OPTIONAL rejected.fasta.gz
# OUTPUT OPTIONAL rejected_1.fasta.gz
# OUTPUT OPTIONAL rejected_2.fasta.gz
# OUTPUT OPTIONAL filter.log
# PARAMETER phred64: "Base quality encoding" TYPE [ n: "Sanger", y: "Illumina v1.3-1.5"] DEFAULT n (Select \"Sanger" for Illumina v1.8+, Sanger, Roche/454, Ion Torrent and PacBio data.)
# PARAMETER input.mode: "Input file format" TYPE [ fq: "FASTQ", fa: "FASTA"] DEFAULT fq (Define the file format of the reads file)
# PARAMETER OPTIONAL min.qual.score: "Minimum quality score" TYPE INTEGER (Filter sequence with at least one quality score below the given value.)
# PARAMETER OPTIONAL min.qual.mean: "Minimum mean quality" TYPE INTEGER (Filter reads with quality score mean below the given value.)
# PARAMETER OPTIONAL min.len: "Minimum length" TYPE INTEGER (Select only reads that are longer than the given value.)
# PARAMETER OPTIONAL max.len: "Maximum length" TYPE INTEGER (Select only reads that are shorter than the given value.)
# PARAMETER OPTIONAL max.gc: "Maximum GC content" TYPE INTEGER (Select only reads that have GC content smaller than the given value.)
# PARAMETER OPTIONAL min.gc: "Minimum GC content" TYPE INTEGER (Select only reads that have GC content larger than the given value.)
# PARAMETER OPTIONAL ns.max.p: "Maximum percentage of Ns" TYPE INTEGER (Filter reads for which the percentage of Ns is higher than the given value.)
# PARAMETER OPTIONAL ns.max.n: "Maximun count of Ns" TYPE INTEGER (Filter reads which have more Ns than the given value.)
# PARAMETER OPTIONAL noniupac: "Remove reads with non-standard characters" TYPE [yes, no] DEFAULT no (Filter out reads with characters other than A, C, G, T or N.)
# PARAMETER OPTIONAL derep: "Type of duplicates to filter" TYPE [0: "none", 1: "exact duplicate", 2: "5-prime duplicate", 3: "3-prime duplicate", 4: "reverse complement exact duplicate", 5:"reverse complement 5-prime/3-prime duplicate"] DEFAULT 0 (Type of duplicates to filter.)
# PARAMETER OPTIONAL derep.min: "Number of allowed duplicates" TYPE INTEGER (This option specifies the number of allowed duplicates. For example, to remove reads that occur more than 5 times, you would specify value 6.)
# PARAMETER OPTIONAL lc.dust: "DUST low complexity threshold" TYPE INTEGER (Use DUST method with the given maximum allowed score, between 0 and 100.)
# PARAMETER OPTIONAL lc.entropy: "ENTROPY low complexity threshold" TYPE INTEGER (Use ENTROPY method with the given minimum entropy value, between 0 and 100.)
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "no", y: "yes"] DEFAULT y (Write a log file)
# PARAMETER OPTIONAL output.mode: "Results to write out" TYPE [ filt: "accepted sequences only", both: "accepted and rejected sequences into separate files"] DEFAULT filt (With this section you can define if the sequences that get filtered out are collected to a separate file) 


# Filter fastq and fasta files based on a number of criteria
# KM, EK, 16-04-2012
# MG, 18-04-2012, added matepair functionality
# KM, 22-05-2012, fixed fastq checking
# AMS 08.01.2014 Paired-end reads now handled by PRINSEQ
# AMS 17.2.2014, gzip outputs

# for later use
# PARAMETER OPTIONAL seq.num: "Maximum number of reads" TYPE INTEGER (Only keep the given number number of reads that pass all other filters.)
# PARAMETER OPTIONAL max.qual.mean: "Maximum mean quality" TYPE INTEGER ( Filter reads with quality score mean above the given value.)
# PARAMETER OPTIONAL max.qual.score: "Maximum quality score" TYPE INTEGER (Filter sequence with at least one quality score above the given value.)

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

# binary
binary.prinseq <- c(file.path(chipster.tools.path, "prinseq", "prinseq-lite.pl" ))

# Parameters
filter.params <- paste("")
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

#if (!is.na(max.qual.score)) {
#	filter.params <- paste(filter.params, "-max_qual_score",  max.qual.score)
#}

if (!is.na(min.qual.mean)) {
	filter.params <- paste(filter.params, "-min_qual_mean",  min.qual.mean)
}

#if (!is.na(max.qual.mean)) {
#	filter.params <- paste(filter.params, "-max_qual_mean",  max.qual.mean)
#}

if (!is.na(ns.max.p)) {
	filter.params <- paste(filter.params, "-ns_max_p",  ns.max.p)
}

if (!is.na(ns.max.n)) {
	filter.params <- paste(filter.params, "-ns_max_n",  ns.max.n)
}

#if (!is.na(seq.num)) {
#	filter.params <- paste(filter.params, "-seq_num",  seq.num)
#}

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

if (input.mode == "fq") {
	filter.params <- paste(filter.params, "-no_qual_header -fastq fastqfile")
	if (is_paired_end) {
		filter.params <- paste(filter.params, "-fastq2 matepair_fastqfile")
	} 
}

if (input.mode == "fa") {
	filter.params <- paste(filter.params, "-fasta fastqfile")
	if (is_paired_end) {
		filter.params <- paste(filter.params, "-fasta2 matepair_fastqfile")
	} 
}

if (output.mode == "both") {
	filter.params <- paste(filter.params, "-out_bad rejected")
}else{
	filter.params <- paste(filter.params, "-out_bad null")
}

filter.command <- paste(binary.prinseq, filter.params, "-out_good accepted")

if (log.file == "y") {
	system("echo Running PRINSEQ filtering with command: >> filter.log")
	echo.command <- paste("echo '", filter.command, "'>> filter.log")
	system(echo.command)
	filter.command <- paste(filter.command, "-verbose 2>> filter.log")
}

system(filter.command)

system("gzip *.fastq")
system("gzip *.fasta")

