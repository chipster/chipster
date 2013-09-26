# TOOL mothur-summaryseqs.R: "Summarize sequences with Mothur" (Given a fasta file with unaligned or aligned sequences, provides summary statistics on sequence start and end coordinates, length, number of ambiguous bases, and homopolymer length. This tool is based on the Mothur package.)
# INPUT reads.fasta: "FASTA file" TYPE FASTA
# OUTPUT summary.tsv
		
# AMS 4.6.2013
# EK 27.6.2013 Changes to description and output
		
# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# batch file
write("summary.seqs(fasta=reads.fasta)", "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log_raw.txt")

# run
system(command)

# Post process output
# system("mv reads.summary summary.tsv")
system("grep -A 9 Start log_raw.txt > summary.tsv")

