# TOOL prinseq-length-filter.R: "Filter reads for length" (Selects nucleotide sequences from a FASTQ file based on read length.)
# INPUT fastqfile: "Input sequence set" TYPE GENERIC
# OUTPUT OPTIONAL accepted.fastq
# OUTPUT OPTIONAL accepted.fasta
# OUTPUT OPTIONAL rejected.fastq
# OUTPUT OPTIONAL rejected.fasta
# OUTPUT OPTIONAL filter.log
# PARAMETER max.len: "Maximum length" TYPE INTEGER (Select only sequences that are shorter than the given value.)
# PARAMETER OPTIONAL min.len: "Minimum length" TYPE INTEGER DEFAULT 0 (Select only sequences that are longer than the given value.)
# PARAMETER OPTIONAL output.mode: "Results to write out" TYPE [ filt: "Accepted sequences only", both: "Accepted and rejected sequences into separate files"] DEFAULT filt (With this section you can define if the sequences that get filtered out are collected to a separate file) 
# PARAMETER OPTIONAL input.mode: "Input file format" TYPE [ fq: "FASTQ", fa: "FASTA"] DEFAULT fq (Define the file format of the reads file)
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "No", y: "Yes"] DEFAULT n (Write a log file. The log file shows the PRINSEQ command used and the amount of sequences in the result files)


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

