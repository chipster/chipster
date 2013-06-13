# TOOL mothur-alignseqs.R: "Align sequences with Mothur" (Align sequences in a fasta-formatted sequence file. This tool is based on the Mothur package.)
# INPUT reads.fasta: "FASTA file" TYPE GENERIC
# OUTPUT batch.mth
# OUTPUT OPTIONAL filtered-aligned.fasta
# OUTPUT OPTIONAL log.tsv


# EK 05.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "1.28.0", "mothur"))

# batch file
write("align.seqs(fasta=reads.fasta, template=/opt/chipster/tools/mothur/data/silva.bacteria.fasta)", "batch.mth", append=F)
write("filter.seqs(fasta=reads.align, vertical=T, trump=.)", "batch.mth", append=T)

# command
command <- paste(binary, "batch.mth")

# run
system(command)

# batch file 2
write("summary.seqs(reads.filter.fasta)", "summary.mth", append=F)

# command
command <- paste(binary, "summary.mth", "> log_raw.txt 2>err.txt")

# run
system(command)

# Post process output
system("grep -A 9 Start log_raw.txt > log.tsv")
system("mv reads.filter.fasta filtered-aligned.fasta")

