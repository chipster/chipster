# TOOL mothur-uniqueseqs.R: "Unique aligned sequences with Mothur" (Unique aligned sequences in a fasta-formatted sequence file. This tool is based on the Mothur package.)
# INPUT a.fasta: "FASTA file" TYPE GENERIC
# INPUT a.names: "Names file" TYPE GENERIC
# OUTPUT OPTIONAL unique.fasta
# OUTPUT OPTIONAL unique.names.txt
# OUTPUT OPTIONAL log.txt


# EK 06.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "1.28.0", "mothur"))

# batch file
write("unique.seqs(fasta=a.fasta, name=a.names)", "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log.txt 2>&1")

system(command)

# Post process output
system("mv a.unique.fasta unique.fasta")
system("mv a.unique.names unique.names.txt")

