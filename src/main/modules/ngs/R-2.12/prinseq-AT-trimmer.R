# TOOL prinseq-AT-trimmer.R: "Trim reads for poly-A/T tails" (Removes poly-A/T tails from reads. This tool is based on the PRINSEQ package.)
# INPUT fastqfile: "Input sequence set" TYPE GENERIC
# OUTPUT OPTIONAL trimmed.fastq.gz
# OUTPUT OPTIONAL trim.log
# PARAMETER OPTIONAL trim.tail.left: "Trim left tails" TYPE INTEGER (Trim poly-A/T tail with a minimum length of the given value at the 5-prime end.)
# PARAMETER OPTIONAL trim.tail.right: "Trim right tails" TYPE INTEGER (Trim poly-A/T tail with a minimum length of the given value at the 3-prime end.)
# PARAMETER OPTIONAL input.mode: "Input file format" TYPE [ fq: "FASTQ", fa: "FASTA"] DEFAULT fq (Define the file format of the reads file)
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "No", y: "Yes"] DEFAULT y (Write a log file)

# KM 17.1.2012
# AMS 17.2.2014, gzip outputs

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("fastqfile")


# binary
binary.prinseq <- c(file.path(chipster.tools.path, "prinseq", "prinseq-lite.pl"))

trim.params <- paste("")

if (!is.na(trim.tail.left)) {
	trim.params <- paste(trim.params, "-trim_tail_left",  trim.tail.left)
}

if (!is.na(trim.tail.right)) {
	trim.params <- paste(trim.params, "-trim_tail_right",  trim.tail.right)
}

if (input.mode == "fq") {
 trim.command <- paste(binary.prinseq, trim.params, "-fastq fastqfile -out_good trimmed -no_qual_header")
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

system("gzip *.fastq")
system("gzip *.fasta")
