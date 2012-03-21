# TOOL fastx-quality-filter.R: "Filter out low quality reads" (Filters out reads which contain quality scores below the user-specified criteria. This tool is based on the FASTQ Quality Filter tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT quality-filtered.fastq
# OUTPUT quality-filtered.log
# PARAMETER quality: "Quality cut-off value" TYPE INTEGER FROM 1 TO 100 DEFAULT 20 (What is the minimum quality score to keep.)
# PARAMETER percentage: "Minimum percent of bases that must have that quality" TYPE INTEGER FROM 1 TO 100 DEFAULT 90 (Percent of bases in sequence that must have quality equal to or higher than the cut-off value.)
# PARAMETER quality.format: "Quality value format used" TYPE [sanger: Sanger, illuminaold: "Illumina GA v1.3-1.5"] DEFAULT sanger (What quality encoding is used in your FASTQ file. Select Sanger if your data comes from Illumina 1.8 or later, SOLiD or 454.)


# EK 28.6.2011

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastq_quality_filter"))

# command
quality.scale <- ifelse(quality.format == "sanger", "-Q 33", "")
command <- paste(binary, "-v", "-q", quality, "-p", percentage, quality.scale, "-i reads.fastq -o quality-filtered.fastq > quality-filtered.log")

# run
system(command)


