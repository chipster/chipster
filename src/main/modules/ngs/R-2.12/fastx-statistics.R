# TOOL fastx-statistics.R: "FASTX quality statistics" (FASTX quality statistics.)
# INPUT reads.fastq TYPE GENERIC 
# OUTPUT quality-stats.tsv 


# EK 17.6.2011

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_quality_stats"))

# command
command <- paste(binary, "-i reads.fastq -o quality-stats.tsv")

# run
system(command)

