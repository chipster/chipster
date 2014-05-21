# TOOL mothur-alignseqs.R: "Align sequences with Mothur" (Given a fasta file of 16S rRNA sequences, aligns them to the Silva reference set. Kmer searching with 8mers is followed by Needleman-Wunsch pairwise alignment, which penalizes the same amount for opening and extending a gap. This tool is based on the Mothur package.)
# INPUT reads.fasta: "FASTA file" TYPE FASTA
# OUTPUT aligned.fasta
# OUTPUT aligned-summary.tsv

# EK 05.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
data.path <- c(file.path(chipster.tools.path, "mothur-data"))
template.path <- c(file.path(data.path, "silva.bacteria.fasta"))

# batch file
write(paste("align.seqs(fasta=reads.fasta, template=", template.path, ")", sep=""), "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth")

# run
system(command)

system("mv reads.align aligned.fasta")

# batch file 2
write("summary.seqs(fasta=aligned.fasta)", "summary.mth", append=F)

# command
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# Post process output
system("grep -A 9 Start log_raw.txt > aligned-summary.tsv")