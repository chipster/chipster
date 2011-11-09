# TOOL fastqc.R: "Read quality with FastQC" (Generates plots for per base quality score and sequence content.
# You can also create plots for GC and N content, sequence length distribution, duplication levels, overrepresented sequences and Kmer content by using the \"Create all plots\" parameter.
# This tool is based on the FastQC package by Simon Andrews et al.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT per_base_quality.png
# OUTPUT per_base_sequence_content.png
# OUTPUT fastqc_data.txt
# OUTPUT OPTIONAL duplication_levels.png
# OUTPUT OPTIONAL kmer_profiles.png
# OUTPUT OPTIONAL per_base_gc_content.png
# OUTPUT OPTIONAL per_base_n_content.png
# OUTPUT OPTIONAL per_sequence_gc_content.png
# OUTPUT OPTIONAL per_sequence_quality.png
# OUTPUT OPTIONAL sequence_length_distribution.png
# PARAMETER all: "Create all plots" TYPE [yes, no] DEFAULT no (Whether to also create plots for GC and N content, sequence length distribution, duplication levels, overrepresented sequences and Kmer content.)

# check out if the file is compressed and if so add a gzip suffix
system("file reads.fastq > file_info")
system("grep gzip file_info > is_gzip")
system("[ -s is_gzip ] && mv reads.fastq reads.gzip")

# binary
binary <- c(file.path(chipster.tools.path, "FastQC", "fastqc"))

# command, different for compressed file
file_info <- scan(file=file_info, nline=1, what="character")
if (grep("gzip", file_info) > 0) {
	command <- paste(binary, "reads.gzip")
} else {
	command <- paste(binary, "reads.fastq")
}

# run
system(command)

# move outputs

system("cp reads_fastqc/Images/per_base_quality.png .")
system("cp reads_fastqc/Images/per_base_sequence_content.png .")
system("cp reads_fastqc/fastqc_data.txt .")

if (all == "yes") {
	system("cp reads_fastqc/Images/duplication_levels.png .")
	system("cp reads_fastqc/Images/kmer_profiles.png .")
	system("cp reads_fastqc/Images/per_base_gc_content.png .")
	system("cp reads_fastqc/Images/per_sequence_quality.png .")
	system("cp reads_fastqc/Images/sequence_length_distribution.png .")
}