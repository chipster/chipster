# TOOL prinseq-N-filter.R: "Filter reads for Ns" (Filters out reads form a reads file based the number or percentage of uniassigned nucleotides in a read.)
# INPUT fastqfile: "Input sequence set" TYPE GENERIC
# OUTPUT OPTIONAL accepted.fastq
# OUTPUT OPTIONAL accepted.fasta
# OUTPUT OPTIONAL rejected.fastq
# OUTPUT OPTIONAL rejected.fasta
# OUTPUT OPTIONAL filter.log
# PARAMETER OPTIONAL ns.max.n: "Maximun count of Ns" TYPE INTEGER (Filter sequences for which the count of Ns id higher than the given value.)
# PARAMETER OPTIONAL ns.max.p: "Maximum percentage of Ns" TYPE INTEGER (Filter sequences for which the percentage of Ns id higher than the given value.)
# PARAMETER OPTIONAL output.mode: "Results to write out" TYPE [ filt: "Accepted sequences only", both: "Accepted and rejected sequences into separate files"] DEFAULT filt (With this section you can define if the sequences that get filtered out are collected to a separate file) 
# PARAMETER OPTIONAL input.mode: "Input file format" TYPE [ fq: "FASTQ", fa: "FASTA"] DEFAULT fq (Define the file format of the reads file)
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "No", y: "Yes"] DEFAULT n (Write a log file)


# check out if the file is compressed and if so unzip it
system("file fastqfile > file_info")
system("grep gzip file_info > is_gzip")
system("[ -s is_gzip ] && mv fastqfile reads.gz ; gzip -d reads.gz ; mv reads fastqfile")


system("
wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.17.3.tar.gz
tar zxf prinseq-lite-0.17.3.tar.gz
")

# binary
#binary.prinseq <- c(file.path(chipster.tools.path, "prinseq", "bin", "prinseq-lite.pl"))
binary.prinseq <- c("perl prinseq-lite-0.17.3/prinseq-lite.pl")



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
