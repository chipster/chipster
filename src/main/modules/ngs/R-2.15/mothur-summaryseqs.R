# TOOL mothur-summaryseqs.R: "Summarize sequences with Mothur" (Summarize the quality of sequences in an unaligned or aligned fasta-formatted sequence file. This tool is based on the Mothur package.)
# INPUT reads.fasta: "FASTA file" TYPE GENERIC
# OUTPUT summary.tsv

# AMS 04.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "1.28.0", "mothur"))

# batch file
write("summary.seqs(fasta=reads.fasta)", "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth")

# run
system(command)

system("mv reads.summary summary.tsv")