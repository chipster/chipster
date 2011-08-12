# TOOL fastqc.R: "FastQC" (FastQC)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT duplication_levels.png
# OUTPUT kmer_profiles.png
# OUTPUT per_base_gc_content.png
# OUTPUT per_base_n_content.png
# OUTPUT per_base_quality.png
# OUTPUT per_base_sequence_content.png
# OUTPUT per_sequence_gc_content.png
# OUTPUT per_sequence_quality.png
# OUTPUT sequence_length_distribution.png

# binary
binary <- c(file.path(chipster.tools.path, "FastQC", "fastqc"))

# command
command <- paste(binary, "reads.fastq")

# run
system(command)

# move outputs
system("cp reads_fastqc/Images/*.png .")

