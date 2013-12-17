# TOOL promoter-retrprom.R: "Retrieve promoters" (Retrieves promoters of the selected genes from UCSC genome database. Currently supports human, mouse, rat, drosophila, and yeast data.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT seqs.txt: seqs.txt 
# PARAMETER species: Species TYPE [human: human, mouse: mouse, rat: rat, drosophila: drosophila, yeast: yeast] DEFAULT human ()
# PARAMETER promoter.size: "Promoter size" TYPE [1000: small, 2000: medium, 5000: large] DEFAULT 1000 (Length of upstream sequences. Small=1000 bp, medium=2000 bp and large=5000 bp)
# PARAMETER multiple.promoters: "Retrieve multiple promoters per gene" TYPE [yes: yes, no: no] DEFAULT no (In the case where the gene has more than one transcription start site, print them all)

# JTT 17.7.2007
# MK 23.09.2013, modified to use R

# Function available in common/R-2.12/promoter.utils. The function returns an object with three vectors (seq, seq.names and nos.names) for sequences, IDs and IDs without sequences
source(file.path(chipster.common.path, "promoter-utils.R"))
seqs <- retreive_promoters(species, promoter.size, multiple.promoters, "normalized.tsv", "phenodata.tsv")

# Retrieving the sequences
for(i in 1:length(seqs$seq.names)) {
      write(file="seqs.txt", seqs$seq.names[i], append=T)
      write(file="seqs.txt", seqs$seq[i], append=T)
}

for(i in 1:length(seqs$nos.names)) {
      write(file="seqs.txt", seqs$nos.names[i], append=T)
      write(file="seqs.txt", paste("", sep=""), append=T)
}
