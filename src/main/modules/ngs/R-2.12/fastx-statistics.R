# TOOL fastx-statistics.R: "Read quality statistics with FASTX" (Calculates quality statistics of the reads. This tool is based on the FASTX Statistics tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC 
# OUTPUT quality-stats.tsv 
# PARAMETER quality.format: "Quality value format used" TYPE [sanger: Sanger, illuminaold: "Illumina GA v1.3-1.5"] DEFAULT sanger (What quality encoding is used in your FASTQ file. Select Sanger if your data comes from Illumina 1.8 or later, SOLiD or 454.)


# EK 17.6.2011

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")

# binary
binary.stats <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_quality_stats"))

# command
quality.scale <- ifelse(quality.format == "sanger", "-Q 33", "")
command.stats <- paste(binary.stats, quality.scale, "-i reads.fastq -o quality-stats.tsv")

# run
ret <- system(command.stats)
if (ret > 0) {
	stop('Unsupported input file type, please see tool output for more details.')
}

# binary
#binary.qualities <- c(file.path(chipster.tools.path, "fastx", "bin", "fastq_quality_boxplot_graph.sh"))

# command
#command.qualities <- paste(binary.qualities, "-i quality-stats.tsv", "-o qualities.png")

# run
#ret <- system(command.qualities)
#if (ret > 0) {
#	stop('Unsupported input file type, please see tool output for more details.')
#}

# binary
#binary.distribution <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_nucleotide_distribution_graph.sh"))

# command
#command.distribution <- paste(binary.distribution, "-i quality-stats.tsv", "-o nucleotide-distribution.png")

# run
#system(command.distribution)

