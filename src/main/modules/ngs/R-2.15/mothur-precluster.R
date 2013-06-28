# TOOL mothur-precluster.R: "Precluster aligned sequences with Mothur" (Preclusters aligned sequences in order to remove sequences that are likely to contain sequencing errors. This tool is based on the Mothur package. In addition to the alignment, you need to supply the names file that was created by the tool \"Extract unique aligned sequences with Mothur\".)
# INPUT a.fasta: "FASTA file" TYPE FASTA
# INPUT a.names: "Names file" TYPE MOTHUR_NAMES
# OUTPUT preclustered.fasta
# OUTPUT preclustered.names
# OUTPUT preclustered-summary.tsv


# EK 18.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# batch file
write("pre.cluster(fasta=a.fasta, name=a.names, diffs=1)", "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log_raw.txt 2>&1")

#run
system(command)

# Post process output
system("mv a.precluster.fasta preclustered.fasta")
system("mv a.precluster.names preclustered.names")

#stool.trim.unique.good.filter.unique.precluster.map

# batch file 2
write("summary.seqs(fasta=preclustered.fasta)", "summary.mth", append=F)

# command 2
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# Post process output
system("grep -A 9 Start log_raw.txt > preclustered-summary.tsv")