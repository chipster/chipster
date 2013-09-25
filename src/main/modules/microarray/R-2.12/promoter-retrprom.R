# TOOL promoter-retrprom.R: "Retrieve promoters" (Retrieves promoters of the selected genes from UCSC genome database. Currently supports human, mouse, rat, drosophila, and yeast data.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT seqs.txt: seqs.txt 
# PARAMETER species: species TYPE [human: human, mouse: mouse, rat: rat, drosophila: drosophila, yeast: yeast] DEFAULT human ()
# PARAMETER promoter.size: promoter.size TYPE [1000: small, 2000: medium, 5000: large] DEFAULT 1000 (Length of upstream sequences. Small=1000 bp, medium=2000 bp and large=5000 bp)
# PARAMETER multiple.promoters: multiple.promoters TYPE [yes: yes, no: no] DEFAULT no (In the case where the gene has more than one transcription start site, print them all)

# JTT 17.7.2007
# MK 23.09.2013, modified to use R

# Load Sequence library
library(GenomicFeatures)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Reads the chiptype from phenodata table
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Loads the correct annotation library
lib<-as.character(phenodata$chiptype[1])
if (length(grep(".db", lib)) == 0 & length(grep("pmcdf", lib)) == 0) {
   lib <- paste(lib, ".db", sep="")
}
library(package=lib, character.only=T)

# Creates a variable for environement
if(species=="drosophila") {
   env<-paste(lib, "ENSEMBL", sep="")
} else {
   env<-paste(lib, "ENTREZID", sep="")
}
env <- sub( ".db", "", env) # if chip contained ".db", remove it

# Retrieve entrezid / ensembl list 
refseq<-as.vector(unlist(mget(rownames(dat), envir=get(env))))
#select(get(lib), rownames(dat), "ENTREZID", "PROBEID")

if(species=="human") {
   genome <-  "BSgenome.Hsapiens.UCSC.hg19"
   database <- "TxDb.Hsapiens.UCSC.hg19.knownGene" 
   organism <- "Hsapiens"
   library(org.Hs.eg.db)
}
if(species=="mouse") {
   genome <-  "BSgenome.Mmusculus.UCSC.mm10"
   database <- "TxDb.Mmusculus.UCSC.mm10.knownGene" 
   organism <- "Mmusculus"
   library(org.Mm.eg.db)
}
if(species=="rat") {
   genome <-  "BSgenome.Rnorvegicus.UCSC.rn5"
   database <- "TxDb.Rnorvegicus.UCSC.rn5.refGene" 
   organism <- "Rnorvegicus"
   library(org.Rn.eg.db)
}
if(species=="drosophila") {
   genome <-  "BSgenome.Dmelanogaster.UCSC.dm3"
   database <- "TxDb.Dmelanogaster.UCSC.dm3.ensGene" 
   organism <- "Dmelanogaster"
   library(org.Dm.eg.db)
}
if(species=="yeast") {
   genome <-  "BSgenome.Scerevisiae.UCSC.sacCer3"
   database <- "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene" 
   organism <- "Scerevisiae"
   library(org.Sc.sgd.db)
}

# Extract Upstream
library(database, character.only=T)
library(genome, character.only=T)

database.entries <- names(transcriptsBy(get(database), by="gene"))
found.refseq <- refseq[!is.na(match(refseq, database.entries))]

if(length(found.refseq) > 0) {
   chromosomal.loc <- transcriptsBy(get(database), by="gene") [found.refseq]
   promoter <- getPromoterSeq(chromosomal.loc, get(organism), upstream=as.integer(promoter.size), downstream=0)
} else {
   stop("CHIPSTER-NOTE: None of your input probes IDs were linked with Refseq IDs")
}


# Retrieving the sequences
if(multiple.promoters=="no") {
   for(i in 1:length(promoter)) {      
      write(file="seqs.txt", paste(">", rownames(dat)[which(refseq==names(promoter)[i])] , " (Refseq ID: ", names(promoter)[i], ")", sep=""), append=T)
      write(file="seqs.txt", paste(toString(promoter[[i]][1]), sep=""), append=T)
   }
} else {
   for(i in 1:length(promoter)) {
      for(j in 1:length(promoter[[i]])) {
         write(file="seqs.txt", paste(">", rownames(dat)[which(refseq==names(promoter)[i])] , " (Refseq ID: ", names(promoter)[i], ")", sep=""), append=T)
         write(file="seqs.txt", paste(toString(promoter[[i]][j]), sep=""), append=T)
      }
   }   
}

#print sequence information for probes without sequence information
for(i in 1:length(rownames(dat))) {
   if(!is.na(match(refseq, database.entries))[i]==FALSE) {
      write(file="seqs.txt", paste(">", rownames(dat)[i], sep=""), append=T)
      write(file="seqs.txt", paste("", sep=""), append=T)
   }
}

