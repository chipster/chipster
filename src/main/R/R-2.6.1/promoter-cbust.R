# ANALYSIS "Promoter Analysis"/ClusterBuster (Does a search for known transcription factor binding sites from
# input sequences using matrices from the Jaspar database. Promoter sequences are automatically retrieved from
# a central database. Currently works only for human, mouse, rat, drosophila, and yeast data.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT clusters.txt
# PARAMETER species [human, mouse, rat, drosophila, yeast] DEFAULT human (Species)
# PARAMETER promoter.size [small, medium, large] DEFAULT small (Length of upstream sequences)
# PARAMETER cluster.score.threshold INTEGER FROM 1 TO 100 DEFAULT 5 (Print details of all motif clusters with score >= this value)
# PARAMETER motif.score.threshold INTEGER FROM 1 TO 100 DEFAULT 6 (Print details of all motif matches that have score >= this value and occur within printed clusters)
# PARAMETER expect.dist.to.neighbor.motifs INTEGER FROM 1 TO 5000 DEFAULT 35 (The expected distance in bp between neighboring motifs in a cluster)
# PARAMETER range.for.counting.nucl.freqs INTEGER FROM 1 TO 10000 DEFAULT 100 (Range in bp for counting local nucleotide abundances)
# PARAMETER pseudocount DECIMAL FROM 0 TO 100 DEFAULT 0.375 (This value gets added to all entries in the motif matrices)


# JTT 15.12.2006

# Sets up the path to the promoter sequences
path.seq<-c("/v/linux26_x86_64/appl/molbio/weeder/seqs/")

# Sets up the path to ClusterBuster executable
path.cbust<-c("/v/linux26_x86_64/appl/molbio/ClusterBuster/cbust-linux")

# Sets up the path to Jaspar Core matrix file
path.jaspar<-c("/v/linux26_x86_64/appl/molbio/ClusterBuster/jaspar2005core.txt")

# Recodes the variable names
size<-promoter.size
cscore<-cluster.score.threshold 
mscore<-motif.score.threshold 
expdist<-expect.dist.to.neighbor.motifs 
nfrange<-range.for.counting.nucl.freqs 
pcount<-pseudocount

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

# Creates a list of genes
genes<-row.names(dat2)

# Loads the annotation library
lib<-as.character(chip)
library(package=lib, character.only=T)

# Creating a list of RefSeq IDs for promoter retrieval
refseq<-as.vector(unlist(mget(genes, envir=get(env))))


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

# Finding known TFBSs from the promoters
arg<-paste("-c ", cscore, " -m ", mscore, " -g ", expdist, " -r ", nfrange, " -p ", pcount, " ", sep="")
system(paste(path.cbust, " ", arg, path.jaspar, " seqs.txt > clusters.txt", sep=""))
