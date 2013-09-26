# TOOL mothur-uniqueseqs.R: "Extract unique aligned sequences with Mothur" (Removes identical sequences from fasta-formatted alignment files. This tool is based on the Mothur package. In addition to the alignment, you need to supply the names file that was created by the tool \"Trim and filter reads with Mothur\". )
# INPUT a.names: "Names file" TYPE MOTHUR_NAMES
# INPUT a.fasta: "FASTA file" TYPE FASTA
# OUTPUT unique.fasta
# OUTPUT unique.names
# OUTPUT unique-summary.tsv


# EK 06.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# batch file
write("unique.seqs(fasta=a.fasta, name=a.names)", "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log.txt 2>&1")

# run
system(command)

# Post process output
system("mv a.unique.fasta unique.fasta")
system("mv a.unique.names unique.names")

# batch file 2
write("summary.seqs(fasta=unique.fasta)", "summary.mth", append=F)

# command 2
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# Post process output
system("grep -A 9 Start log_raw.txt > unique-summary.tsv")