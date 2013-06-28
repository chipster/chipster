# TOOL mothur-classifyseqs.R: "Classify sequences to taxonomic units with Mothur" (Classify 16S rRNA sequences to taxonomic units using the Wang method with Silva reference set and taxonomy. This tool is based on the Mothur package.)
# INPUT a.fasta: "FASTA file" TYPE FASTA
# OUTPUT reads-taxonomy-assignment.txt
# OUTPUT classification-summary.txt


# EK 18.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
data.path <- c(file.path(chipster.tools.path, "mothur-data"))
template.path <- c(file.path(data.path, "silva.bacteria.fasta"))
taxonomy.path <- c(file.path(data.path, "silva.bacteria.silva.tax"))

# batch file
write(paste("classify.seqs(fasta=a.fasta, iters=1000, template=", template.path, ", taxonomy=", taxonomy.path, ")", sep=""), "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log.txt 2>&1")

# run
system(command)

# Post process output
system("mv a.silva.wang.taxonomy reads-taxonomy-assignment.txt")
system("mv a.silva.wang.tax.summary classification-summary.txt")

