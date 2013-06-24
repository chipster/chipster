# TOOL mothur-filterseqs.R: "Filter aligned sequences with Mothur" (Filter aligned sequences in a fasta-formatted sequence file. This tool is based on the Mothur package.)
# INPUT a.align: "Aligned reads in FASTA format" TYPE GENERIC
# OUTPUT filtered-aligned.fasta
# OUTPUT OPTIONAL log.txt

# EK 05.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "1.28.0", "mothur"))

# batch file
write("filter.seqs(fasta=a.align, vertical=T, trump=.)", "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log_raw.txt")

# run
system(command)

system("mv a.filter.fasta filtered-aligned.fasta")

system("grep -A 4 Length log_raw.txt > log.txt")