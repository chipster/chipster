# TOOL convert-bam-to-edger.R: "Convert BAM file to edgeR input format" (This tool takes BAM files as an input, calculates the number of times each sequence tag is identified and removes the ones for which the count is under the user defined threshold)
# INPUT bam_file.bam: "BAM data dile" TYPE GENERIC
# OUTPUT edgeR-input{...}.tsv: "A converted BAM file suitable for edgeR analysis"
# PARAMETER count_limit: "Count limit" TYPE INTEGER FROM 0 TO 1000 DEFAULT 10 (The lowest number of times a sequence tag has to appear in the data)

# MG 15.6.2011

# Sam tools setup
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
command.start <- paste("bash -c '", samtools.binary)

# Extract the BAM file into SAM
samtools.parameters <- "view"
samtools.input <- "bam_file.bam"
samtools.output <- "sam_file.view"
samtools.command <- paste(command.start, " ", samtools.parameters, " ", samtools.input,
				" > ", samtools.output)
system(samtools.command)

samtools bam_file.bam > sam_file_view

# Extract the following info from the SAM file:
# nucleotide sequence
# chromosome
# start
# length of sequence
awk '{print $10"\t"$3"\t"$4"\t"length($10)+$4}' miRNA.sorted_view > miRNA.sorted_col_3_4_10_TEST


# Find the unique reads
uniq -c -f 1 miRNA.sorted_col_3_4_10_TEST > miRNA.sorted_unique


# sort them
sort -k3n -k4n miRNA.sorted_unique > miRNA.sorted_sorted

# make a usable output file
awk '{if($1> 10)print $2"\t"$3"\t"$4"\t"$1}' miRNA.sorted_sorted > peaksHigherThan10

# remove sterisk chromosome
grep -v \* peaksHigherThan10 > peaksHigherThan10_trimmed




# common parameters
common.parameters <- paste("-q --best -S --strata", "-m", multiread, "-k", alignment.no)

# mode specific parameters
quality.parameter <- ifelse(quality.format == "solexa1_3", "--solexa1.3-quals", "")
n.mode.parameters <- paste("-n", max.mismatches, "-l", seed, "-e", quality, quality.parameter)
v.mode.parameters <- paste("-v", max.mismatches)
mode.parameters <- ifelse(limit.to.seed == "yes", n.mode.parameters, v.mode.parameters)

# output parameters
unaligned.output <- ifelse(unaligned.file == "yes", "--un unaligned-reads.fastq", "")
multiread.output <- ifelse(multiread.file == "yes", "--max multireads.fastq", "")
output.parameters <- paste(unaligned.output, multiread.output)

# command ending
command.end <- paste(genome, "reads.txt 1> alignment.sam 2> bowtie.log'")

# run bowtie
bowtie.command <- paste(command.start, common.parameters, mode.parameters, output.parameters, command.end)
#stop(paste('CHIPSTER-NOTE: ', bowtie.command))
system(bowtie.command)


# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# convert sam to bam
system(paste(samtools.binary, "view -bS alignment.sam -o alignment.bam"))

# sort bam
system(paste(samtools.binary, "sort alignment.bam alignment.sorted"))

# index bam
system(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
system("mv alignment.sorted.bam bowtie.bam")
system("mv alignment.sorted.bam.bai bowtie.bam.bai")
