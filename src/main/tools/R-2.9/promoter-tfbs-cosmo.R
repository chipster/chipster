# ANALYSIS "Promoter Analysis"/Cosmo (Finds common sequence motifs in the promoters of input genes. 
# Promoter sequences are automatically retrieved from a central database. Currently works only for human, mouse, 
# rat, drosophila, and yeast data.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT cosmo-output.txt, seqlogo.png, probs.png
# PARAMETER species [human, mouse, rat, drosophila, yeast] DEFAULT human (Species)
# PARAMETER promoter.size [small, medium, large] DEFAULT small (Length of upstream sequences)
# PARAMETER strands [single, both] DEFAULT single (Analyze both strands of DNA)
# PARAMETER appears.more.than.once [yes, no] DEFAULT no (Could the motif appear more than once in every sequence)
# PARAMETER percentage INTEGER FROM 1 TO 100 DEFAULT 50 (Percentage of sequences the motif should appear)
# PARAMETER tfsize [small, medium] DEFAULT small (Transcription factor binding site size)


# Promoter sequence analysis
# JTT 21.5.2008

# Loads the libraries
library(cosmo)

# Sets up the path to the promoter sequences
path.seq<-c(file.path(chipster.tools.path, "weeder", "seqs"))

# Renaming variable
size<-promoter.size
once<-appears.more.than.once

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Read phenodata and extracts chip information
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
chip<-phenodata$chiptype[1]

# Creates a variable for environement
if(species=="drosophila") {
   env<-paste(chip, "ACCNUM", sep="")
} else {
   env<-paste(chip, "REFSEQ", sep="")
}
env <- sub( ".db", "", env) # if chip contained ".db", remove it

# Creates a list of genes
genes<-row.names(dat)

# Loads the annotation library
lib<-as.character(chip)

# Account for the fact that annotation packages are from version 2.3 of Bioconductor
# named with an ".db" suffix. Add the suffix when missing to support data files
# from Chipster 1.3 and earlier. 
if (length(grep(".db", lib)) == 0 & length(grep("pmcdf", lib)) == 0) {
        lib <- paste(lib, ".db", sep="")
}

library(package=lib, character.only=T)

# Creating a list of RefSeq IDs for promoter retrieval
refseq<-as.vector(unlist(mget(genes, envir=get(env))))
refseq<-unique(refseq)

# Retrieving promoters
if(species=="human" & size=="small") {
   upstream<-read.table(paste(path.seq, "UCSC_Build_36.1_hg18_upstream1000.tsv", sep=""), header=T, sep="\t")
}
if(species=="human" & size=="medium") {
   upstream<-read.table(paste(path.seq, "UCSC_Build_36.1_hg18_upstream2000.tsv", sep=""), header=T, sep="\t")
}
if(species=="human" & size=="large") {
   upstream<-read.table(paste(path.seq, "UCSC_Build_36.1_hg18_upstream5000.tsv", sep=""), header=T, sep="\t")
}
if(species=="mouse" & size=="small") {
   upstream<-read.table(paste(path.seq, "UCSC_Build_36_mm8_upstream1000.tsv", sep=""), header=T, sep="\t")
}
if(species=="mouse" & size=="medium") {
   upstream<-read.table(paste(path.seq, "UCSC_Build_36_mm8_upstream2000.tsv", sep=""), header=T, sep="\t")
}
if(species=="mouse" & size=="large") {
   upstream<-read.table(paste(path.seq, "UCSC_Build_36_mm8_upstream5000.tsv", sep=""), header=T, sep="\t")
}
if(species=="rat" & size=="small") {
   upstream<-read.table(paste(path.seq, "UCSC_rn4_upstream1000.tsv", sep=""), header=T, sep="\t")
}
if(species=="rat" & size=="medium") {
   upstream<-read.table(paste(path.seq, "UCSC_rn4_upstream2000.tsv", sep=""), header=T, sep="\t")
}
if(species=="rat" & size=="large") {
   upstream<-read.table(paste(path.seq, "UCSC_rn4_upstream5000.tsv", sep=""), header=T, sep="\t")
}
if(species=="drosophila" & size=="small") {
   upstream<-read.table(paste(path.seq, "Drosophila_upstream1000.tsv", sep=""), header=T, sep="\t")
}
if(species=="drosophila" & size=="medium") {
   upstream<-read.table(paste(path.seq, "Drosophila_upstream2000.tsv", sep=""), header=T, sep="\t")
}
if(species=="drosophila" & size=="large") {
   upstream<-read.table(paste(path.seq, "Drosophila_upstream5000.tsv", sep=""), header=T, sep="\t")
}
if(species=="yeast" & size=="small") {
   upstream<-read.table(paste(path.seq, "NCBI_sc_upstream500.tsv", sep=""), header=T, sep="\t")
}
if(species=="yeast" & size=="medium") {
   upstream<-read.table(paste(path.seq, "NCBI_sc_upstream1000.tsv", sep=""), header=T, sep="\t")
}
if(species=="yeast" & size=="large") {
   upstream<-read.table(paste(path.seq, "NCBI_sc_upstream2500.tsv", sep=""), header=T, sep="\t")
}

# Retrieving the sequences
w<-c()
for(i in 1:length(refseq)) {
   w<-c(w, which(upstream$RefSeq==refseq[i]))
}
unlink("seqs.txt")
for(i in 1:length(w)) {
   write(file="seqs.txt", paste(">", upstream[w[i],]$RefSeq, sep=""), append=T)
   write(file="seqs.txt", paste(upstream[w[i],]$Sequence, sep=""), append=T)
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

bitmap(file="seqlogo.png", width=600/72, height=600/72)
plot(res)
dev.off()

bitmap(file="probs.png", width=600/72, height=600/72)
plot(res, type="prob")
dev.off()