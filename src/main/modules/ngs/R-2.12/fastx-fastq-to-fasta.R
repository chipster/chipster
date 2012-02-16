# TOOL fastx-fastq-to-fasta.R: "Convert FASTQ to FASTA" (Convert FASTQ files to FASTA format. This tool is based on the FASTX package.)
# INPUT reads.fastq TYPE GENERIC 
# OUTPUT reads.fasta
# OUTPUT fasta.log
# PARAMETER OPTIONAL remove.unknowns: "Remove sequences with unknown nucleotides" TYPE [yes, no] DEFAULT yes (Remove sequences with unknown nucleotides)
# PARAMETER OPTIONAL rename.identifiers: "Rename sequence identifiers as numbers" TYPE [yes, no] DEFAULT no (Rename sequence identifiers as numbers)



# EK 17.6.2011

# check out if the file is compressed and if so unzip it
system("file reads.fastq > file_info")
system("grep gzip file_info > is_gzip")
system("[ -s is_gzip ] && mv reads.fastq reads.gz ; gzip -d reads.gz ; mv reads reads.fastq")

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastq_to_fasta"))

# parameters
remove.parameter <- ifelse(remove.unknowns == "yes", "", "-n")
rename.parameter <- ifelse(rename.identifiers == "no", "", "-r")

# command
command <- paste(binary, remove.parameter, rename.parameter, "-v -i reads.fastq -o reads.fasta > fasta.log")

# run
system(command)

