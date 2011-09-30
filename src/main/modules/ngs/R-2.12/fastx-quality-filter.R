# TOOL fastx-quality-filter.R: "Filter out low quality reads" (Filters out reads which contain quality scores below the user-specified criteria. This tool is based on the FASTQ Quality Filter tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT quality-filtered.fastq
# OUTPUT quality-filtered.log
# PARAMETER quality: "Quality cut-off value" TYPE INTEGER FROM 1 TO 100 DEFAULT 5 (What is the minimum quality score to keep.)
# PARAMETER percentage: "Minimum percent of bases that must have that quality" TYPE INTEGER FROM 1 TO 100 DEFAULT 90 (Percent of bases in sequence that must have quality equal to or higher than the cut-off value.)



# EK 28.6.2011

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastq_quality_filter"))

# command
command <- paste(binary, "-v", "-q", quality, "-p", percentage, "-i reads.fastq -o quality-filtered.fastq > quality-filtered.log")

# run
system(command)


