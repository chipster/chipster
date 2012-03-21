# TOOL prinseq-AT-trimmer.R: "Poly-A/T trimmer" (Removes poly-A/T tails from nucleotide sequences. This tool utilizes PRINSEQ program.)
# INPUT fastqfile: "Input sequence set" TYPE GENERIC
# OUTPUT trimmed.fastq
# OUTPUT OPTIONAL trim.log
# PARAMETER OPTIONAL trim.tail.left: "Trim left tails" TYPE INTEGER (Trim poly-A/T tail with a minimum length of the given value at the 5-prime-end.)
# PARAMETER OPTIONAL trim.tail.right: "Trim right tails" TYPE INTEGER (Trim poly-A/T tail with a minimum length of the given value at the 3-prime-end.)
# PARAMETER OPTIONAL input.mode: "Input file format" TYPE [ fq: "FASTQ", fa: "FASTA"] DEFAULT fq (Define the file format of the reads file)
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "No", y: "Yes"] DEFAULT n (Write a log file)



# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("fastqfile")

system("
wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.17.3.tar.gz
tar zxf prinseq-lite-0.17.3.tar.gz
")

# binary
#binary.prinseq <- c(file.path(chipster.tools.path, "prinseq", "bin", "prinseq-lite.pl"))
binary.prinseq <- c("perl prinseq-lite-0.17.3/prinseq-lite.pl")



trim.params <- paste(" ")

if (!is.na(trim.tail.left)) {
	trim.params <- paste(trim.params, "-trim_tail_left",  trim.tail.left)
}

if (!is.na(trim.tail.right)) {
	trim.params <- paste(trim.params, "-trim_tail_right",  trim.tail.right)
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
	system("if [ ! -e  trimmed.fasta ] ; then echo 'Trimming produced an empty trimmed.fasta sequence set' >> trim.log ; echo '' > trimmed.fasta ; fi")
}


#stop(paste('CHIPSTER-NOTE: ', filter.command))

