# TOOL mothur-chimera-uchime.R: "Remove chimeric sequences with Mothur" (Remove chimeric sequences from a fasta-formatted file. This tool is based on the Mothur package.)
# INPUT a.fasta: "FASTA file" TYPE GENERIC
# OUTPUT OPTIONAL chimeras-removed.fasta
# OUTPUT OPTIONAL log.txt


# EK 18.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "1.28.0", "mothur"))
data.path <- c(file.path(chipster.tools.path, "mothur", "data"))
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

system("grep -A 2 Removed log_raw.txt > log.txt")

