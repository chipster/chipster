# TOOL fastx-collapser.R: "Remove duplicate reads from FASTQ or FASTA" (Collapses identical sequences in a FASTQ/FASTA file into a single sequence. The sequences are renamed with sequence number and the multiplicity value. This tool is based on the FASTA/Q Collapser tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT duplicates-removed.fastq 


# EK 12.1.2012

# check out if the file is compressed and if so unzip it
system("file reads.fastq > file_info")
system("grep gzip file_info > is_gzip")
system("[ -s is_gzip ] && mv reads.fastq reads.gz ; gzip -d reads.gz ; mv reads reads.fastq")

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_collapser"))

# command
command <- paste(binary, "-i reads.fastq -o duplicates-removed.fastq")

# run
system(command)