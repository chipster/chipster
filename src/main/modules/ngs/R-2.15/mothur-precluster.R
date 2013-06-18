# TOOL mothur-precluster.R: "Precluster aligned sequences with Mothur" (Precluster aligned sequences in a fasta-formatted sequence file. This tool is based on the Mothur package.)
# INPUT a.fasta: "FASTA file" TYPE GENERIC
# INPUT a.names: "Names file" TYPE GENERIC
# OUTPUT OPTIONAL preclustered.fasta
# OUTPUT OPTIONAL preclustered.names.txt
# OUTPUT OPTIONAL log.txt


# EK 18.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "1.28.0", "mothur"))

# batch file
write("pre.cluster(fasta=a.fasta, name=a.names, diffs=1)", "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log_raw.txt 2>&1")

#run
system(command)

# Post process output
system("mv a.precluster.fasta preclustered.fasta")
system("mv a.precluster.names preclustered.names.txt")

system("grep -A 2 Total log_raw.txt > log.txt")

#stool.trim.unique.good.filter.unique.precluster.map