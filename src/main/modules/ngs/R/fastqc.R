# TOOL fastqc.R: "Read quality with FastQC" (Generates plots for per base quality score and sequence content, read GC and N content, length distribution, duplication levels, overrepresented sequences and Kmer content.
# This tool is based on the FastQC package by Simon Andrews et al.)
# INPUT reads TYPE GENERIC
# OUTPUT fastqc_report.pdf
# PARAMETER filetype: "File type" TYPE [fastq: "FASTQ", bam: "BAM"] DEFAULT fastq (Select input file type.)

# 2014.12.16 AMS Changed output to PDF, removed parameter for all plots

library(png)
library(gplots)

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


# Create pdf
#
pdf("fastqc_report.pdf", width=210/25.4, height=297/25.4) # page size A4
par(mar=c(0,0,3,0), mfrow=c(2, 1))

#page 1
summary <- read.table("reads_fastqc/summary.txt", header=FALSE, , sep = "\t")
system("csplit -s reads_fastqc/fastqc_data.txt /END_MODULE/ \"{*}\"") #split the result file into separate parts for each analysis module (xx00..xx11)
basics <- read.table("xx00", header=FALSE, skip=3, sep = "\t")
textplot(summary, show.rownames=FALSE, show.colnames=FALSE, halign="center", valign="top", cex=0.80)
title("Summary")
textplot(basics, show.rownames=FALSE, show.colnames=FALSE, halign="center", valign="top", cex=0.80)
title("Basic statistics")

#page 2
plot.new()
img <- readPNG("reads_fastqc/Images/per_base_quality.png", native = TRUE, info = FALSE)
rasterImage(img, 0, 0, 1, 1, interpolation=TRUE)
title("Per base sequence quality")
plot.new()
img <- readPNG("reads_fastqc/Images/per_sequence_quality.png", native = TRUE, info = FALSE)
rasterImage(img, 0, 0, 1, 1, interpolation=TRUE)
title("Per sequence quality scores")

#page 3
plot.new()
img <- readPNG("reads_fastqc/Images/per_base_sequence_content.png", native = TRUE, info = FALSE)
rasterImage(img, 0, 0, 1, 1, interpolation=TRUE)
title("Per base sequence content")
plot.new()
img <- readPNG("reads_fastqc/Images/per_base_gc_content.png", native = TRUE, info = FALSE)
rasterImage(img, 0, 0, 1, 1, interpolation=TRUE)
title("Per base GC content")

#page 4
plot.new()
img <- readPNG("reads_fastqc/Images/per_sequence_gc_content.png", native = TRUE, info = FALSE)
rasterImage(img, 0, 0, 1, 1, interpolation=TRUE)
title("Per sequence GC content")
plot.new()
img <- readPNG("reads_fastqc/Images/per_base_n_content.png", native = TRUE, info = FALSE)
rasterImage(img, 0, 0, 1, 1, interpolation=TRUE)
title("Per base N content")

#page 5
plot.new()
img <- readPNG("reads_fastqc/Images/sequence_length_distribution.png", native = TRUE, info = FALSE)
rasterImage(img, 0, 0, 1, 1, interpolation=TRUE)
title("Sequence Length Distribution")
plot.new()
img <- readPNG("reads_fastqc/Images/duplication_levels.png", native = TRUE, info = FALSE)
rasterImage(img, 0, 0, 1, 1, interpolation=TRUE)
title("Sequence Duplication Levels")

#page 6
if (file.info("xx09")$size > 50){
	system("sed s/#// xx09 > ors.txt") # get rid of hash sign in col name
	ors <- read.table("ors.txt", header=TRUE, skip=2, sep = "\t")
	textplot(ors, show.rownames=FALSE, show.colnames=TRUE, halign="left", valign="top", cex=0.7)
	title("Overrepresented sequences")
}

#page 7
plot.new()
img <- readPNG("reads_fastqc/Images/kmer_profiles.png", native = TRUE, info = FALSE)
rasterImage(img, 0, 0, 1, 1, interpolation=TRUE)
title("Kmer Content")
dev.off()



