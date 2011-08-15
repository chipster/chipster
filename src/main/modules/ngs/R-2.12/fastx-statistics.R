# TOOL fastx-statistics.R: "Read quality statistics with FASTX" (Calculates quality statistics of the reads. Please note that this tool only works with Sanger quality scores.
# This tool is based on the FASTX Statistics tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC 
# OUTPUT quality-stats.tsv 


# EK 17.6.2011

# binary
binary.stats <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_quality_stats"))

# command
command.stats <- paste(binary.stats, "-i reads.fastq -o quality-stats.tsv")

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

