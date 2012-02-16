# TOOL fastx-trimmer.R: "Trim reads" (Trims reads to a user-specified length. This tool is based on the FASTA/Q Trimmer tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT trimmed.fastq 
# PARAMETER first: "First base to keep" TYPE INTEGER FROM 1 TO 100 DEFAULT 1 (First base to keep.)
# PARAMETER last: "Last base to keep" TYPE INTEGER FROM 1 TO 300 DEFAULT 75 (Last base to keep.)




# EK 17.6.2011

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_trimmer"))

# check out if the file is compressed and if so unzip it
system("file reads.fastq > file_info")
system("grep gzip file_info > is_gzip")
system("[ -s is_gzip ] && mv reads.fastq reads.gz ; gzip -d reads.gz ; mv reads reads.fastq")

# command
command <- paste(binary, "-f", first, "-l", last, "-i reads.fastq -o trimmed.fastq")

# run
system(command)