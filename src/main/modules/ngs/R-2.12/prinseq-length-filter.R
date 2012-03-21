# TOOL prinseq-length-filter.R: "Filter reads based on length" (Filters reads based on length. This tool is based on the PRINSEQ package.)
# INPUT fastqfile: "Input sequence set" TYPE GENERIC
# OUTPUT OPTIONAL accepted.fastq
# OUTPUT OPTIONAL accepted.fasta
# OUTPUT OPTIONAL rejected.fastq
# OUTPUT OPTIONAL rejected.fasta
# OUTPUT OPTIONAL filter.log
# PARAMETER max.len: "Maximum length" TYPE INTEGER DEFAULT 500 (Select only reads that are shorter than the given value.)
# PARAMETER min.len: "Minimum length" TYPE INTEGER DEFAULT 15 (Select only reads that are longer than the given value.)
# PARAMETER OPTIONAL output.mode: "Results to write out" TYPE [ filt: "accepted reads only", both: "accepted and rejected reads into separate files"] DEFAULT filt (With this section you can define if the reads that get filtered out are collected to a separate file) 
# PARAMETER OPTIONAL input.mode: "Input file format" TYPE [ fq: "FASTQ", fa: "FASTA"] DEFAULT fq (Define the file format of the reads file)
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "no", y: "yes"] DEFAULT n (Write a log file. The log file shows the PRINSEQ command used and the amount of reads in the result files)


# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("fastqfile")

# binary
binary.prinseq <- c(file.path(chipster.tools.path, "prinseq", "prinseq-lite.pl"))



filter.params <- paste(" ")
filter.params <- paste(filter.params, "-max_len",  max.len)

if (!is.na(min.len)) {
   filter.params <- paste(filter.params, "-min_len",  min.len)
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
	system("if [ ! -e  accepted.fastq ] ; then echo 'Filtering produced an empty accepted.fastq sequence set.' >> filter.log ; echo '' > accepted.fastq ; fi")
}

if (input.mode == "fa") {
	system("if [ ! -e  accepted.fasta ] ; then echo 'Filtering produced an empty accepted.fasta sequence set.' >> filter.log  ; echo '' > accepted.fasta ; fi")
}

if (output.mode == "both") {
	if (input.mode == "fq") {
		system("if [ ! -e  rejected.fastq ] ; then echo 'Filtering produced an empty rejected.fastq sequence set.' >> filter.log ; echo '' > rejected.fastq ; fi")
	}
	
	if (input.mode == "fa") {
		system("if [ ! -e  rejected.fasta ] ; then echo 'Filtering produced an empty rejected.fasta sequence set.' >> filter.log  ; echo '' > rejected.fasta ; fi")
	}
}


#stop

