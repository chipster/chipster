# ANALYSIS "Promoter Analysis"/"Cosmo with Ensembl" (Finds common sequence motifs in the promoters of input genes. 
# Promoter sequences are automatically retrieved from a central database. Currently works only for human and mouse for ensembl and 
# rat, drosophila, and yeast data as well for ucsc. If using Ensembl, use promotor.size for cosmo data and upstream/downstream for ensembl promotors. Modified by Mandy.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT cosmo-output.txt, seqlogo.png, probs.png
# PARAMETER species [human, mouse, rat, drosophila, yeast] DEFAULT human (Species)
# PARAMETER promoter.size [small, medium, large] DEFAULT small (Length of upstream sequences)
# PARAMETER strands [single, both] DEFAULT single (Analyze both strands of DNA)
# PARAMETER appears.more.than.once [yes, no] DEFAULT no (Could the motif appear more than once in every sequence)
# PARAMETER percentage INTEGER FROM 1 TO 100 DEFAULT 50 (Percentage of sequences the motif should appear)
# PARAMETER tfsize [small, medium] DEFAULT small (Transcription factor binding site size)
# PARAMETER Database.Selection [ucsc, ensembl] DEFAULT ucsc (Database to query)
# PARAMETER upstreamc INTEGER FROM 1 TO 2000 DEFAULT 500 (How many upstream)
# PARAMETER downstreamc INTEGER FROM 0 TO 1500 DEFAULT 0 (How many downstream)


# Promoter sequence analysis
# JTT 21.5.2008
# MbyM 16.2.09

# Loads the libraries
library(cosmo)

# Sets up the path to the promoter sequences
path.seq<-c("/home/s2686739/Desktop/Dataforchipster/promodata/")

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
if (Database.Selection=="ensembl"){
	env<-chip
}else{
	if(species=="drosophila") {
  		 env<-paste(chip, "ACCNUM", sep="")
	} else {
   		env<-paste(chip, "REFSEQ", sep="")
	}
}

# Creates a list of genes
genes<-row.names(dat)

# Loads the annotation library
lib<-as.character(chip)
library(package=lib, character.only=T)

if (Database.Selection=="ucsc"){
	refseq<-as.vector(unlist(mget(genes, envir=get(env))))
	refseq<-unique(refseq)
}else{
	refseq<-row.names(dat)
}
# Getting data from ensembl files
if(Database.Selection=="ensembl"){
	if (species=="mouse"){
		upstream<-read.table(paste(path.seq, "Ensemblpromodata/MouseEnsembl2000up1700down.tsv", sep=""), header=T, sep="\t")	
	}else{
		upstream<-read.table(paste(path.seq, "Ensemblpromodata/HumanEnsembl2000up1700down.tsv", sep=""), header=T, sep="\t")
	}
	w<-c()
	for(i in 1:length(refseq)) {
		if (species=="human"){
   			w<-c(w, which(upstream@row.names==refseq[i]))
		}else{
			w<-c(w, which(upstream[,1]==refseq[i]))
		}
	}
	#making seq file
	for(i in 1:length(w)) {
		up<-(2000 - upstreamc)
		down<-(2000 + downstreamc)
		print(w)
		#unlink("seqs.txt")
		if (species=="mouse"){
			d<-(upstream[w[i],2])
  			write(file="seqs.txt", paste(">", upstream[w[i],1], sep=""), append=TRUE)
   			write(file="seqs.txt", paste(substr(d,up,down), sep=""), append=TRUE)

		}else{
			d<-(upstream[w[i],1])
  		 	write(file="seqs.txt", paste(">", upstream@row.names[w[i]], sep=""), append=TRUE)
   			write(file="seqs.txt", paste(substr(d,up,down), sep=""), append=TRUE)
		}
	}

}else{

	# Retrieving promoters
	if(species=="human" & size=="small") {
	   upstream<-read.table(paste(path.seq, "hg18/UCSC_Build_36.1_hg18_upstream1000.tsv", sep=""), header=T, sep="\t")
	}
	if(species=="human" & size=="medium") {
	   upstream<-read.table(paste(path.seq, "hg18/UCSC_Build_36.1_hg18_upstream2000.tsv", sep=""), header=T, sep="\t")
	}
	if(species=="human" & size=="large") {
	   upstream<-read.table(paste(path.seq, "hg18/UCSC_Build_36.1_hg18_upstream5000.tsv", sep=""), header=T, sep="\t")
	}
	if(species=="mouse" & size=="small") {
	   upstream<-read.table(paste(path.seq, "mm8/UCSC_Build_36_mm8_upstream1000.tsv", sep=""), header=T, sep="\t")
	}
	if(species=="mouse" & size=="medium") {
	   upstream<-read.table(paste(path.seq, "mm8/UCSC_Build_36_mm8_upstream2000.tsv", sep=""), header=T, sep="\t")
	}
	if(species=="mouse" & size=="large") {
	   upstream<-read.table(paste(path.seq, "mm8/UCSC_Build_36_mm8_upstream5000.tsv", sep=""), header=T, sep="\t")
	}
	if(species=="rat" & size=="small") {
	   upstream<-read.table(paste(path.seq, "Rat/UCSC_rn4_upstream1000.tsv", sep=""), header=T, sep="\t")
	}
	if(species=="rat" & size=="medium") {
	   upstream<-read.table(paste(path.seq, "Rat/UCSC_rn4_upstream2000.tsv", sep=""), header=T, sep="\t")
	}
	if(species=="rat" & size=="large") {
	   upstream<-read.table(paste(path.seq, "Rat/UCSC_rn4_upstream5000.tsv", sep=""), header=T, sep="\t")
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
	   upstream<-read.table(paste(path.seq, "yeast/NCBI_sc_upstream500.tsv", sep=""), header=T, sep="\t")
	}
	if(species=="yeast" & size=="medium") {
	   upstream<-read.table(paste(path.seq, "yeast/NCBI_sc_upstream1000.tsv", sep=""), header=T, sep="\t")
	}
	if(species=="yeast" & size=="large") {
	   upstream<-read.table(paste(path.seq, "yeast/NCBI_sc_upstream2500.tsv", sep=""), header=T, sep="\t")
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
