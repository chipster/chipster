# TOOL convert-bam-to-edger.R: "Convert BAM file to edgeR input format" (This tool takes BAM files as an input, calculates the number of times each sequence tag is identified and removes the ones for which the count is under the user defined threshold)
# INPUT bam_file.bam: "BAM data dile" TYPE GENERIC
# OUTPUT edgeR-input.tsv: "A converted BAM file suitable for edgeR analysis"
# PARAMETER count_limit: "Count limit" TYPE INTEGER FROM 0 TO 1000 DEFAULT 10 (The lowest number of times a sequence tag has to appear in the data)

# MG 15.6.2011

# Sam tools setup
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
#command.start <- paste("bash -c '", samtools.binary)
#command.start <- "/fs/local/users/chipster/tools/samtools/samtools"
command.start <- samtools.binary
# Extract the BAM file into SAM
samtools.parameters <- "view"
samtools.input <- "bam_file.bam"
samtools.output <- "sam_file"
samtools.command <- paste(command.start, samtools.parameters, samtools.input,
				">", samtools.output)
system(samtools.command)

# Extract the following info from the SAM file:
# nucleotide sequence
# chromosome
# start
# length of sequence
input.file <- "sam_file"
output.file <- "sam_file_extracted"
extract.command <- paste ("awk '{print $10\"\t\"$3\"\t\"$4\"\t\"length($10)+$4-1\"\t\"length($10)}'", input.file, ">", output.file)
system(extract.command)


# Sort sequence reads according to chromosome and start position
input.file <- "sam_file_extracted"
output.file <- "sam_file_sorted"
sort.command <- paste ("sort -k3n -k4n", input.file, ">", output.file)
system(sort.command)

# Find the unique reads
input.file <- "sam_file_sorted"
output.file <- "sam_file_unique"
unique.command <- paste ("uniq -c -f 1", input.file, ">", output.file)
system(unique.command)

# Create an output file with sequence reads that occur at least count_limit times
input.file <- "sam_file_unique"
output.file <- "sam_file_output"
output.command <- paste ("awk '{if($1>", count_limit, ")print $2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$1}'", input.file, ">", output.file)
system(output.command)

# Creat sequence read ID composed of chromosome name, start position and end position
input.file <- "sam_file_output"
output.file <- "sam_file_id"
id.command <- paste ("awk '{print $2\"_\"$3\"_\"$4\"_\"$1\"\t\"$1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}'", input.file, ">", output.file)
system(id.command)

# Remove sequence reads mapping to random chromosome
input.file <- "sam_file_id"
output.file <- "sam_file_trimmed"
trim.command <- paste ("grep -v \\*", input.file, ">", output.file)
system(trim.command)

# Add column headers
headers <- paste("id\t","sequence\t","chr\t","start\t","end\t","length\t","count", sep="")
input.file <- "sam_file_trimmed"
header.file <- "header_file"
output.file <- "edgeR-input.tsv"
system(paste("echo \"", headers, "\"", ">", header.file))
merge.command <- paste("cat", header.file, input.file, ">", output.file)
system(merge.command)

# EOF

