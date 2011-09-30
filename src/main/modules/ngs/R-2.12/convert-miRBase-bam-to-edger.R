# TOOL convert-miRBase-bam-to-edger.R: "Convert miRBased BAM file to edgeR input format" (This tool takes BAM files as an input, calculates the number of times each sequence tag is identified and removes the ones for which the count is under the user defined threshold.)
# INPUT bam_file.bam: "Alignment agains miRBase in BAM format" TYPE GENERIC
# OUTPUT edgeR-input.tsv: "A converted BAM file suitable for edgeR analysis"
# PARAMETER count_limit: "Count limit" TYPE INTEGER FROM 0 TO 1000 DEFAULT 10 (Keep miRNAs which have more reads than this.)

# EK 7.7.2011

# Convert to SAM, grep miRNA name, count, filter on tag number
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
counts.command <- paste(samtools.binary, "view bam_file.bam | awk '{print $3}' | grep -v '^*' | sort -k1 | uniq -c | awk '{if($1>", count_limit, ")print $2\"\t\"$1}'> counts.tsv")
system(counts.command)


# Add column headers
headers <- paste("id\t","count", sep="")
input.file <- "counts.tsv"
header.file <- "header_file"
system(paste("echo \"", headers, "\"", ">", header.file))
merge.command <- paste("cat", header.file, input.file, "> edgeR-input.tsv")
system(merge.command)

#headers <- paste("id\t","count", sep="")
#system(paste("echo \"", headers, "\"", "> header.file"))
#merge.command <- paste("cat", header.file, "counts.tsv > edgeR-input.tsv")
#system(merge.command)

# EOF

