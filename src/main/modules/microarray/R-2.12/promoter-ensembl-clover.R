# ANALYSIS "Promoter Analysis"/"Clover with Ensembl" (Clover Cis eLement OVERrepresentation, see cagt.bu.edu/page/Clover_about for more information. Includes Ensembl database for Human and Mouse. Made by Mandy.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT outputclover.txt
# PARAMETER species [human, mouse, rat, drosophila, yeast] DEFAULT human (Species)
# PARAMETER size [small, medium, large] DEFAULT small (Length of upstream sequences)
# PARAMETER RCrawS INTEGER FROM 0 TO 100 DEFAULT 0 (Number of randomized/control raw scores to calculate for comparison with each target raw score. Setting zero will disclude this function)
# PARAMETER Pvalue INTEGER FROM 0 TO 100 DEFAULT 0 (P-value threshold: only print results for motifs whose P-values don't exceed this amount. Setting zero will disclude this function)
# PARAMETER STPrintingLoc INTEGER FROM 0 TO 200 DEFAULT 0 (Score threshold for printing locations of significant motifs. This parameter doesn't affect raw score and P-value calculations, just which motif instances get printed. Setting zero will disclude this function)
# PARAMETER nucShuffles [yes, no] DEFAULT no (Perform sequence shuffles)
# PARAMETER DenucRandomization [yes, no] DEFAULT no (Perform dinucleotide randomizations.)
# PARAMETER MotifShuffles [yes, no] DEFAULT no (Perform Motif Shuffles.)
# PARAMETER MaskLowercase [yes, no] DEFAULT no (Mask any lowercase letters in the target sequences and background sequences, if any. Lowercase letters are often used to indicate repetitive elements)
# PARAMETER Database.Selection [ucsc] DEFAULT ucsc (Database to query)
# PARAMETER upstreamc INTEGER FROM 1 TO 2000 DEFAULT 500 (How many upstream)
# PARAMETER downstreamc INTEGER FROM 0 TO 1500 DEFAULT 0 (How many downstream)
# PARAMETER Motifs [jaspar] DEFAULT jaspar ( What form you would like the motifs to be in, Jasper MAO or regular)

# MbyM 2-2-09

# Sets up the path to the promoter sequences
path.seq<-c("/home/s2686739/Desktop/Dataforchipster/promodata/")

# Sets up the path to ClusterBuster executable
path.cbust<-c("/home/s2686739/Desktop/Dataforchipster/clover")


# Sets up the path to Jaspar Core matrix file
path.jaspar<-c("/home/s2686739/Desktop/Dataforchipster/jaspar2005core.txt")
if(Motifs=="regular"){
	path.jaspar<-c("/home/s2686739/Desktop/Dataforchipster/jaspar2005coredifnames.txt")
}


# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Read phenodata and extracts chip information
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
chip<-phenodata$chiptype[1]
if (Database.Selection=="ensembl"){
	env<-chip
}else{
	if(species=="drosophila") {
  		 env<-paste(chip, "ACCNUM", sep="")
	} else {
   		env<-paste(chip, "REFSEQ", sep="")
	}
}
env <- sub( ".db", "", env) # if chip contained ".db", remove it

genes<-row.names(dat2)
lib<-as.character(chip)

# Account for the fact that annotation packages are from version 2.3 of Bioconductor
# named with an ".db" suffix. Add the suffix when missing to support data files
# from Chipster 1.3 and earlier. 
if (length(grep(".db", lib)) == 0 & length(grep("pmcdf", lib)) == 0) {
        lib <- paste(lib, ".db", sep="")
}

library(package=lib, character.only=T)
if (Database.Selection=="ucsc"){
	refseq<-as.vector(unlist(mget(genes, envir=get(env))))
	refseq<-unique(refseq)
}else{
	refseq<-row.names(dat)
}
if(Database.Selection=="ensembl"){
	if (species=="mouse"){
		upstream<-read.table(file.path(path.seq, "Ensemblpromodata/MouseEnsembl2000up1700down.tsv"), header=T, sep="\t")	
	}else{
		upstream<-read.table(file.path(path.seq, "Ensemblpromodata/HumanEnsembl2000up1700down.tsv"), header=T, sep="\t")
	}
	w<-c()
	for(i in 1:length(refseq)) {
		if (species=="human"){
   			w<-c(w, which(upstream@row.names==refseq[i]))
		}else{
			w<-c(w, which(upstream[,1]==refseq[i]))
		}
	}
	# Making Seq files
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

	if(species=="human" & size=="small") {
   		upstream<-read.table(file.path(path.seq, "hg18/UCSC_Build_36.1_hg18_upstream1000.tsv"), header=T, sep="\t")
	}
	if(species=="human" & size=="medium") {
   		upstream<-read.table(file.path(path.seq, "hg18/UCSC_Build_36.1_hg18_upstream2000.tsv"), header=T, sep="\t")
	}
	if(species=="human" & size=="large") {
   		upstream<-read.table(file.path(path.seq, "hg18/UCSC_Build_36.1_hg18_upstream5000.tsv"), header=T, sep="\t")
	}
	if(species=="mouse" & size=="small") {
   		upstream<-read.table(file.path(path.seq, "mm8/UCSC_Build_36_mm8_upstream1000.tsv"), header=T, sep="\t")
	}
	if(species=="mouse" & size=="medium") {
  		 upstream<-read.table(file.path(path.seq, "mm8/UCSC_Build_36_mm8_upstream2000.tsv"), header=T, sep="\t")
	}
	if(species=="mouse" & size=="large") {
   		upstream<-read.table(file.path(path.seq, "mm8/UCSC_Build_36_mm8_upstream5000.tsv"), header=T, sep="\t")
	}
	if(species=="rat" & size=="small") {
   		upstream<-read.table(file.path(path.seq, "Rat/UCSC_rn4_upstream1000.tsv"), header=T, sep="\t")
	}
	if(species=="rat" & size=="medium") {
   		upstream<-read.table(file.path(path.seq, "Rat/UCSC_rn4_upstream2000.tsv"), header=T, sep="\t")
	}
	if(species=="rat" & size=="large") {
   		upstream<-read.table(file.path(path.seq, "Rat/UCSC_rn4_upstream5000.tsv"), header=T, sep="\t")
	}
	if(species=="drosophila" & size=="small") {
   		upstream<-read.table(file.path(path.seq, "Drosophila_upstream1000.tsv"), header=T, sep="\t")
	}
	if(species=="drosophila" & size=="medium") {
   		upstream<-read.table(file.path(path.seq, "Drosophila_upstream2000.tsv"), header=T, sep="\t")
	}
	if(species=="drosophila" & size=="large") {
   		upstream<-read.table(file.path(path.seq, "Drosophila_upstream5000.tsv"), header=T, sep="\t")
	}
	if(species=="yeast" & size=="small") {
   		upstream<-read.table(file.path(path.seq, "yeast/NCBI_sc_upstream500.tsv"), header=T, sep="\t")
	}
	if(species=="yeast" & size=="medium") {
   		upstream<-read.table(file.path(path.seq, "yeast/NCBI_sc_upstream1000.tsv"), header=T, sep="\t")
	}
	if(species=="yeast" & size=="large") {
   		upstream<-read.table(file.path(path.seq, "yeast/NCBI_sc_upstream2500.tsv"), header=T, sep="\t")
	}
	w<-c()
	for(i in 1:length(refseq)) {
  		w<-c(w, which(upstream$RefSeq==refseq[i]))
	
	}
	#unlink("seqs.txt")
	for(i in 1:length(w)) {
  		write(file="seqs.txt", paste(">", upstream[w[i],]$RefSeq, sep=""), append=T)
   		write(file="seqs.txt", paste(upstream[w[i],]$Sequence, sep=""), append=T)
		print(upstream[w[i],]$RefSeq)
	}
}
args<-" "
if (nucShuffles=="yes"){
	args<-paste(args," -n ")
	}
if (DenucRandomization=="yes"){
	args<-paste(args," -d ")
	}
if (MotifShuffles=="yes"){
	args<-paste(args," -m ")
	}
if (MaskLowercase=="yes"){
	args<-paste(args," -l ")
	}
if (RCrawS>0){
args<-paste(args," -r ", RCrawS," ",sep="")
}
if (Pvalue>0){
args<-paste(args," -t ", Pvalue," ",sep="")
}
if (STPrintingLoc>0){
args<-paste(args," -u ", STPrintingLoc," ",sep="")
}
system("head -n 20 seqs.txt")
system(paste(path.cbust, args, path.jaspar, " seqs.txt > outputclover.txt", sep=""))

