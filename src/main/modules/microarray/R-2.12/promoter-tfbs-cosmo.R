# TOOL promoter-tfbs-cosmo.R: Cosmo (Finds common sequence motifs in the promoters of input genes. Promoter sequences are automatically retrieved from a central database. Currently works only for human, mouse, rat, drosophila, and yeast data.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT cosmo-output.txt: cosmo-output.txt 
# OUTPUT seqlogo.pdf: seqlogo.pdf 
# OUTPUT probs.pdf: probs.pdf 
# PARAMETER species: Species TYPE [human: human, mouse: mouse, rat: rat, drosophila: drosophila, yeast: yeast] DEFAULT human ()
# PARAMETER promoter.size: "Promoter size" TYPE [1000: small, 2000: medium, 5000: large] DEFAULT 1000 (Length of upstream sequences. Small=1000 bp, medium=2000 bp and large=5000 bp)
# PARAMETER multiple.promoters: "Retrieve multiple promoters per gene" TYPE [yes: yes, no: no] DEFAULT no (In the case where the gene has more than one transcription start site, print them all)
# PARAMETER strands: Strands TYPE [single: single, both: both] DEFAULT single (Analyze both strands of DNA)
# PARAMETER appears.more.than.once: "Appears more than once" TYPE [yes: yes, no: no] DEFAULT no (Could the motif appear more than once in every sequence)
# PARAMETER percentage: Percentage TYPE INTEGER FROM 1 TO 100 DEFAULT 50 (Percentage of sequences the motif should appear)
# PARAMETER tfsize: tfsize TYPE [small: small, medium: medium] DEFAULT small (Transcription factor binding site size)


# JTT 21.05.2008
# MK 25.09.2013: promoter sequence retrieve modified
# MK 01.10.2013: script removed from the tool category, since R.3.0.0 does not support cosmo anymore


# Loads the libraries
library(cosmo)

# Sets up the path to the promoter sequences
path.seq<-c(file.path(chipster.tools.path, "weeder", "seqs"))

# Renaming variable
size<-promoter.size
once<-appears.more.than.once

# Retrieving the sequences. The function is available in common/R-2.12/promoter.utils.
source(file.path(chipster.common.path, "promoter-utils.R"))
seqs <- retreive_promoters(species, promoter.size, multiple.promoters, "normalized.tsv", "phenodata.tsv")

# Write sequences on disk
for(i in 1:length(seqs$seq.names)) {
      write(file="seqs.txt", seqs$seq.names[i], append=T)
      write(file="seqs.txt", seqs$seq[i], append=T)
}

# Generating the function call
cll<-paste("cosmo(seqs='seqs.txt', ", sep="")

if(strands=="single") {
   cll<-paste(cll, "revComp=F, ", sep="")
} else {
   cll<-paste(cll, "revComp=T, ", sep="")
}

if(once=="yes") {
   cll<-paste(cll, "models=c('TCM'), ", sep="")
} else {
   cll<-paste(cll, "models=c('ZOOPS'), ", sep="")
}

pers<-floor(length(w)*percentage/100)
cll<-paste(cll, "minSites=", pers, ", ", sep="")

if(tfsize=="small") {
   cll<-paste(cll, "minW=6, maxW=8", sep="")
} else {
   cll<-paste(cll, "minW=10, maxW=10", sep="")
}

cll<-paste(cll, ")", sep="")
cll<-paste("res<-", cll, sep="")

# Running the analysis
write(cll, "cosmo-batch.txt")
source("cosmo-batch.txt") 

# Saving the results
sink("cosmo-output.txt")
summary(res)
sink()

pdf(file="seqlogo.pdf", width=600/72, height=600/72)
plot(res)
dev.off()

pdf(file="probs.pdf", width=600/72, height=600/72)
plot(res, type="prob")
dev.off()
