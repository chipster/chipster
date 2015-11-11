# TOOL fastqc.R: "Read quality with FastQC" (Generates plots for per base quality score and sequence content, read GC and N content, length distribution, duplication levels, overrepresented sequences and Kmer content.
# This tool is based on the FastQC package by Simon Andrews et al.)
# INPUT reads TYPE GENERIC
# OUTPUT reads_fastqc.html
# PARAMETER filetype: "File type" TYPE [fastq: "FASTQ", bam: "BAM"] DEFAULT fastq (Select input file type.)

# 2014.12.16 AMS Changed output to PDF, removed parameter for all plots
# 2015.09.10 AMS New version embeds pictures in html, so changed output to html 

#library(png)
#library(gplots)

# FastQC detects gzipped files by file extension so we need to add .gz
# extension to compressed files.
source(file.path(chipster.common.path, "zip-utils.R"))
input.file <- "reads"
if (isGZipFile(input.file)) {
	system(paste("mv", input.file, "reads.gz"))
	input.file <- "reads.gz"
}

# binary
binary <- file.path(chipster.tools.path, "FastQC", "fastqc")

# command
command <- paste(binary, "-f", filetype, input.file, "--extract")

# run
system(command)

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))


# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("reads_fastqc.html", paste(strip_name(inputnames$reads), "_fastqc.html", sep=""))

# Write output definitions file
write_output_definitions(outputnames)
