# TOOL mothur-alignseqs.R: "Align sequences with Mothur" (Align sequences in a fasta-formatted sequence file. This tool is based on the Mothur package.)
# INPUT reads.fasta: "FASTA file" TYPE GENERIC
# OUTPUT aligned.fasta
# OUTPUT OPTIONAL log.tsv

# EK 05.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "1.28.0", "mothur"))

# batch file
write("align.seqs(fasta=reads.fasta, template=/opt/chipster/tools/mothur/data/silva.bacteria.fasta)", "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth")

# run
system(command)

system("mv reads.align aligned.fasta")

# batch file 2
# write("summary.seqs(reads.filter.fasta)", "summary.mth", append=F)

# command
# command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
# system(command2)

# Post process output
# system("grep -A 9 Start log_raw.txt > log.tsv")

