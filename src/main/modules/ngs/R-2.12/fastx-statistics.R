# TOOL fastx-statistics.R: "FASTX quality statistics" (FASTX quality statistics.)
# INPUT reads.fastq TYPE GENERIC 
# OUTPUT quality-stats.tsv 
# OUTPUT qualities.png



# EK 17.6.2011

# binary
binary.stats <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_quality_stats"))

# command
command.stats <- paste(binary.stats, "-i reads.fastq -o quality-stats.tsv")

# run
system(command.stats)

# binary
binary.qualities <- c(file.path(chipster.tools.path, "fastx", "bin", "fastq_quality_boxplot_graph.sh"))

# command
command.qualities <- paste(binary.qualities, "-i quality-stats.tsv", "-o qualities.png")

# run
system(command.qualities)

# binary
#binary.distribution <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_nucleotide_distribution_graph.sh"))

# command
#command.distribution <- paste(binary.distribution, "-i quality-stats.tsv", "-o nucleotide-distribution.png")

# run
#system(command.distribution)

