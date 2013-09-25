# TOOL mothur-chimera-uchime.R: "Remove chimeric sequences with Mothur" (Remove chimeric sequences from a fasta-formatted alignment using the uchime method and the 16S rRNA Silva gold reference set. This tool is based on the Mothur package.)
# INPUT a.fasta: "FASTA file" TYPE FASTA
# OUTPUT OPTIONAL chimeras-removed.fasta
# OUTPUT OPTIONAL chimeras-removed-summary.tsv


# EK 18.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
data.path <- c(file.path(chipster.tools.path, "mothur-data"))
template.path <- c(file.path(data.path, "silva.gold.align"))

# batch file
write(paste("chimera.uchime(fasta=a.fasta, template=", template.path, ")", sep=""), "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log_raw.txt 2>&1")

# run
system(command)

# batch file 2
write("remove.seqs(accnos=a.uchime.accnos, fasta=a.fasta)", "remove.mth", append=F)

# command
command2 <- paste(binary, "remove.mth", "> log_raw.txt 2>&1")

# run
system(command2)

# Post process output
system("mv a.pick.fasta chimeras-removed.fasta")

# system("grep -A 2 Removed log_raw.txt > log.txt")

# batch file 3
write("summary.seqs(fasta=chimeras-removed.fasta)", "summary.mth", append=F)

# command 3
command3 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command3)

# Post process output
system("grep -A 9 Start log_raw.txt > chimeras-removed-summary.tsv")