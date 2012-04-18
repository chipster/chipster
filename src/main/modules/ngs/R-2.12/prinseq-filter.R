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
# PARAMETER OPTIONAL phred64: "Quality data is in Phred+64 format" TYPE [ n: "no", y: "yes"] DEFAULT n (You should select \"yes\" option if the quality data in FASTQ file is in Phred+64 format. For Illumina 1.8+, Sanger, Roche/454, Ion Torrent, PacBio data, you should use the default value: no)
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "no", y: "yes"] DEFAULT y (Write a log file)

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("fastqfile")


# Check if two files were given as input and if so run the python script
# that interlaces the mate pairs into a single file
input_files <- dir()
is_paired_end <- (length(grep("matepair_fastqfile", input_files))>0)
if (is_paired_end) {
	# binary
#	binary_python_scripts <- file.path(chipster.module.path, "shell", "match-mate-pairs")
	binary_python_scripts <- file.path("/opt/chipster4/comp/modules/ngs/shell/match-mate-pairs", "interleave-fastq.py")
	system_command <- paste("python", binary_python_scripts, "fastqfile", "matepair_fastqfile", "interleaved_fastqfile")
	system(system_command)	
	system("echo Executed interleave python script with: > filter.log")
	echo.command <- paste("echo '", system_command, "'>> filter.log")
	system(echo.command)
	system("ls -l >> filter.log")
	# remove input files to clear up disk space
	system("rm -f fastqfile")
	system("rm -f matepair_fastqfile")
	system("mv interleaved_fastqfile fastqfile")
	system("ls -l >> filter.log")
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

system("ls -l >> filter.log")
system("echo Beginning of file: >> filter.log")
system("head accepted.fastq >> filter.log")
system("echo Emd of file >> filter.log")
system("tail accepted.fastq >> filter.log")

# If filtering on paired-end data perform matching of
# nate pairs using python script and then de-interlace
if (is_paired_end) {
#	system_command <- paste(binary_python_scripts, "match_pairs.py", "accepted.fastq", "matched_fastqfile")
	binary_python_scripts <- file.path("/opt/chipster4/comp/modules/ngs/shell/match-mate-pairs", "match-pairs.py")
	system_command <- paste("python", binary_python_scripts, "accepted.fastq", "matched_fastqfile")
	system(system_command)
	
	system("echo Executed match_pair python script with: >> filter.log")
	echo.command <- paste("echo '", system_command, "'>> filter.log")
	system(echo.command)

	# remove input files to clear up disk space
	system("rm -f accepted.fastq")
	
	system("ls -l >> filter.log")

#	system_command <- paste(binary_python_scripts, "deinterleave.py", "matched_fastqfile", "accepted.fastq", "accepted_matepair.fastq")
	binary_python_scripts <- file.path("/opt/chipster4/comp/modules/ngs/shell/match-mate-pairs", "deinterleave-fastq.py")
	system_command <- paste("python", binary_python_scripts, "matched_fastqfile", "accepted.fastq", "accepted_matepair.fastq")
	system(system_command)	

	system("echo Executed deinterleave python script with: >> filter.log")
	echo.command <- paste("echo '", system_command, "'>> filter.log")
	system(echo.command)

	system("ls -l >> filter.log")
	
	system("rm -f matched_fastqfile")
}

# stop
