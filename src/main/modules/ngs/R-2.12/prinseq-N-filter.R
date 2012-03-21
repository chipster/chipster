# TOOL prinseq-N-filter.R: "Filter reads based on Ns" (Filters out reads form a reads file based the number or percentage of unassigned nucleotides in a read. This tool is based on the PRINSEQ package.)
# INPUT fastqfile: "Input sequence set" TYPE GENERIC
# OUTPUT OPTIONAL accepted.fastq
# OUTPUT OPTIONAL accepted.fasta
# OUTPUT OPTIONAL rejected.fastq
# OUTPUT OPTIONAL rejected.fasta
# OUTPUT OPTIONAL filter.log
# PARAMETER OPTIONAL ns.max.n: "Maximum count of Ns" TYPE INTEGER (Filter out reads for which the count of Ns is higher than the given value.)
# PARAMETER OPTIONAL ns.max.p: "Maximum percentage of Ns" TYPE INTEGER (Filter reads for which the percentage of Ns id higher than the given value.)
# PARAMETER OPTIONAL output.mode: "Results to write out" TYPE [ filt: "accepted reads only", both: "accepted and rejected reads into separate files"] DEFAULT filt (With this section you can define if the reads that get filtered out are collected to a separate file) 
# PARAMETER OPTIONAL input.mode: "Input file format" TYPE [ fq: "FASTQ", fa: "FASTA"] DEFAULT fq (Define the file format of the reads file)
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "no", y: "yes"] DEFAULT n (Write a log file)


# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("fastqfile")


# binary
binary.prinseq <- c(file.path(chipster.tools.path, "prinseq", "prinseq-lite.pl" ))

filter.params <- paste(" ")

if (!is.na(ns.max.p)) {
	filter.params <- paste(filter.params, "-ns_max_p",  ns.max.p)
}

if (!is.na(ns.max.n))  {
	filter.params <- paste(filter.params, "-ns_max_n",  ns.max.n)
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
	system("echo Running PRINSEQ filtering with command: > filter.log")
	echo.command <- paste("echo '", filter.command, "'>> filter.log")
	system(echo.command)
	filter.command <- paste(filter.command, "-verbose 2>> filter.log")
}

system(filter.command)

#Make sure that someting is given as an output
if (input.mode == "fq") {
	system("if [ ! -e  accepted.fastq ] ; then echo 'Filtering produced an empty accepted.fastq sequence set' >> filter.log ; echo '' > accepted.fastq ; fi")
}

if (input.mode == "fa") {
	system("if [ ! -e  accepted.fasta ] ; then echo 'Filtering produced an empty accepted.fasta sequence set' >> filter.log  ; echo '' > accepted.fasta ; fi")
}

if (output.mode == "both") {
	if (input.mode == "fq") {
		system("if [ ! -e  rejected.fastq ] ; then echo 'Filtering produced an empty rejected.fastq sequence set' >> filter.log ; echo '' > rejected.fastq ; fi")
	}
	
	if (input.mode == "fa") {
		system("if [ ! -e  rejected.fasta ] ; then echo 'Filtering produced an empty rejected.fasta sequence set' >> filter.log  ; echo '' > rejected.fasta ; fi")
	}
}

#stop
