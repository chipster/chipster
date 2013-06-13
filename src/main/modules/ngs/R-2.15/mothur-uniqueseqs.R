# TOOL mothur-uniqueseqs.R: "Unique aligned sequences with Mothur" (Unique aligned sequences in a fasta-formatted sequence file. This tool is based on the Mothur package.)
# INPUT a.fasta: "FASTA file" TYPE GENERIC
# INPUT a.names: "Names file" TYPE GENERIC
# OUTPUT OPTIONAL unique.fasta
# OUTPUT OPTIONAL unique.names.txt


# EK 06.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "1.28.0", "mothur"))

# batch file
write("unique.seqs(fasta=a.fasta, names=a.names)", "batch.mth", append=F)

# command and run
command <- paste(binary, "batch.mth")
system(command)

# Post process output
system("mv a.unique.fasta unique.fasta")
system("mv a.unique.names unique.names.txt")

