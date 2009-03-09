# ANALYSIS "Promoter Analysis"/"Retrieve promoters using Ensembl genes" (Retrieves promoters of the selected genes from Ensembl 
# database.Currently works only for human and mouse. Made by Mandy)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT seqs.txt
# PARAMETER species [human, mouse] DEFAULT human (Species)
# PARAMETER upstreamc INTEGER FROM 1 TO 2000 DEFAULT 500 (How many upstream)
# PARAMETER downstreamc INTEGER FROM 0 TO 1500 DEFAULT 0 (How many downstream)

# JTT 17.11.2006
# MbyM 11-02-09

# Sets up the path to the promoter sequences

path.seq<-c("/home/s2686739/Desktop/Dataforchipster/promodata/Ensemblpromodata/")


# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Read phenodata and extracts chip information
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
chip<-phenodata$chiptype[1]
env<-chip

# Creates a list of genes
genes<-row.names(dat)

# Loads the annotation library
lib<-as.character(chip)
library(package=lib, character.only=T)
refseq<-row.names(dat)
if (species=="mouse"){
	upstream<-read.table(paste(path.seq, "MouseEnsembl2000up1700down.tsv", sep=""), header=T, sep="\t")	
}else{
	upstream<-read.table(paste(path.seq, "HumanEnsembl2000up1700down.tsv", sep=""), header=T, sep="\t")

}

# Retrieving the sequences
w<-c()
for(i in 1:length(refseq)) {
	if (species=="human"){
   		w<-c(w, which(upstream@row.names==refseq[i]))
	}else{
		w<-c(w, which(upstream[,1]==refseq[i]))
	}
}
	#making seq file
unlink("seqs.txt")
for(i in 1:length(w)) {

	up<-(2000 - upstreamc)
	down<-(2000 + downstreamc)
	if (species=="mouse"){
		d<-(upstream[w[i],2])
  		write(file="seqs.txt", paste(">", upstream[w[i],1], sep=""), append=T)
   		write(file="seqs.txt", paste(substr(d,up,down), sep=""), append=T)
	}else{
		d<-(upstream[w[i],1])
  		 write(file="seqs.txt", paste(">", upstream@row.names[w[i]], sep=""), append=T)
   		write(file="seqs.txt", paste(substr(d,up,down), sep=""), append=T)
	}
}


