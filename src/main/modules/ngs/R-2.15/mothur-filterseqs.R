# TOOL mothur-filterseqs.R: "Filter sequence alignment with Mothur" (Filters out columns from a fasta formatted sequence alignment. By removing empty columns, the distance calculation is accelerated. This tool is based on the Mothur package.)
# INPUT a.align: "Aligned reads in FASTA format" TYPE FASTA
# OUTPUT filtered-aligned.fasta
# OUTPUT filtered-log.txt
# OUTPUT filtered-summary.tsv

# EK 05.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# batch file
write("filter.seqs(fasta=a.align, vertical=T, trump=.)", "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log_raw.txt")

# run
system(command)

# result post-processing
system("mv a.filter.fasta filtered-aligned.fasta")
system("grep -A 4 filtered log_raw.txt > filtered-log.txt")

# batch file 2
write("summary.seqs(fasta=filtered-aligned.fasta)", "summary.mth", append=F)

# command 2
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# Post process output
system("grep -A 9 Start log_raw.txt > filtered-summary.tsv")