# ANALYSIS "Promoter Analysis"/Weeder (Finds common sequence motifs in the promoters of input genes. 
# Promoter sequences are automatically retrieved from a central database. Currently works only for human, mouse, 
# rat, drosophila, and yeast data.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT seqs.html, seqs.txt.wee
# PARAMETER species [human, mouse, rat, drosophila, yeast] DEFAULT human (Species)
# PARAMETER promoter.size [small, medium, large] DEFAULT small (Length of upstream sequences)
# PARAMETER strands [single, both] DEFAULT single (Analyze both strands of DNA)
# PARAMETER appears.more.than.once [yes, no] DEFAULT no (Could the motif appear more than once in every sequence)
# PARAMETER no.of.motifs INTEGER FROM 1 TO 100 DEFAULT 10 (Number of motifs to return)
# PARAMETER percentage INTEGER FROM 1 TO 100 DEFAULT 50 (Percentage of sequences the motif should appear)
# PARAMETER tfsize [small, medium] DEFAULT small (Transcription factor binding site size)

# Promoter sequence analysis
# JTT 17.11.2006


# Sets up the path to the promoter sequences
path.seq<-c(file.path(chipster.tools.path, "weeder", "seqs"))

# Sets up the paths to weeder executables
path.weeder.tfbs<-c(file.path(chipster.tools.path, "weeder", "Weeder1.3", "weederTFBS.out"))
path.weeder.advicer<-c(file.path(chipster.tools.path, "weeder", "Weeder1.3", "adviser.out"))

# Sets up the paths to weeder frequency files
path.weeder.freq<-c(file.path(chipster.tools.path, "weeder", "Weeder1.3", "FreqFiles"))

# Renaming variable
size<-promoter.size
once<-appears.more.than.once
no<-no.of.motifs

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

# Weeder needs frequency files - copying those in place
system("mkdir FreqFiles")
if(species=="human") {
   system(paste("ln -s ", path.weeder.freq, "HS.6.freq FreqFiles/HS.6.freq", sep=""))
   system(paste("ln -s ", path.weeder.freq, "HS.8.freq FreqFiles/HS.8.freq", sep=""))
   system(paste("ln -s ", path.weeder.freq, "HSI.6.freq FreqFiles/HSI.6.freq", sep=""))
   system(paste("ln -s ", path.weeder.freq, "HSI.8.freq FreqFiles/HSI.8.freq", sep=""))
}
if(species=="mouse") {
   system(paste("ln -s ", path.weeder.freq, "MM.6.freq FreqFiles/MM.6.freq", sep=""))
   system(paste("ln -s ", path.weeder.freq, "MM.8.freq FreqFiles/MM.8.freq", sep=""))
   system(paste("ln -s ", path.weeder.freq, "MMI.6.freq FreqFiles/MMI.6.freq", sep=""))
   system(paste("ln -s ", path.weeder.freq, "MMI.8.freq FreqFiles/MMI.8.freq", sep=""))
}
if(species=="rat") {
   system(paste("ln -s ", path.weeder.freq, "RN.6.freq FreqFiles/RN.6.freq", sep=""))
   system(paste("ln -s ", path.weeder.freq, "RN.8.freq FreqFiles/RN.8.freq", sep=""))
}
if(species=="drosophila") {
   system(paste("ln -s ", path.weeder.freq, "DM.6.freq FreqFiles/RN.6.freq", sep=""))
   system(paste("ln -s ", path.weeder.freq, "DM.8.freq FreqFiles/RN.8.freq", sep=""))
}
if(species=="yeast") {
   system(paste("ln -s ", path.weeder.freq, "SC.6.freq FreqFiles/SC.6.freq", sep=""))
   system(paste("ln -s ", path.weeder.freq, "SC.8.freq FreqFiles/SC.8.freq", sep=""))
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

# Assembling Weeder arguments
arg<-c()
if(percentage!=50) {
   arg<-paste(arg, " -R ", percentage, sep="")
} else {
   arg<-paste(arg, " -R ", 50, sep="")
}
if(strands=="both") {
   arg<-paste(arg, " -S ", sep="")
}
if(once=="yes") {
   arg<-paste(arg, "-M ", sep="")
}
if(no!=10) {
   arg<-paste(arg, "-T ", no, sep="")
}


# Running Weeder to find tentative TFBSs
if(species=="human" & tfsize=="small") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O HS -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O HS -W 8 -e 2", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}
if(species=="human" & tfsize=="medium") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O HS -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O HS -W 8 -e 2", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O HS -W 10 -e 3", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}
if(species=="human" & tfsize=="large") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O HS -W 12 -e 4", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}
if(species=="human" & tfsize=="extra") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O HS -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O HS -W 8 -e 3", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O HS -W 10 -e 4", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O HS -W 12 -e 4", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}


if(species=="mouse" & tfsize=="small") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O MM -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O MM -W 8 -e 2", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}
if(species=="mouse" & tfsize=="medium") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O MM -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O MM -W 8 -e 2", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O MM -W 10 -e 3", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}
if(species=="mouse" & tfsize=="large") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O MM -W 12 -e 4", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}
if(species=="mouse" & tfsize=="extra") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O MM -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O MM -W 8 -e 3", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O MM -W 10 -e 4", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O MM -W 12 -e 4", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}


if(species=="rat" & tfsize=="small") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O RN -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O RN -W 8 -e 2", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}
if(species=="rat" & tfsize=="medium") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O RN -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O RN -W 8 -e 2", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O RN -W 10 -e 3", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}
if(species=="rat" & tfsize=="large") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O RN -W 12 -e 4", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}
if(species=="rat" & tfsize=="extra") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O RN -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O RN -W 8 -e 3", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O RN -W 10 -e 4", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O RN -W 12 -e 4", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}


if(species=="yeast" & tfsize=="small") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O SC -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O SC -W 8 -e 2", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}
if(species=="yeast" & tfsize=="medium") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O SC -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O SC -W 8 -e 2", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O SC -W 10 -e 3", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}
if(species=="yeast" & tfsize=="large") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O SC -W 12 -e 4", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}
if(species=="yeast" & tfsize=="extra") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O SC -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O SC -W 8 -e 3", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O SC -W 10 -e 4", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O SC -W 12 -e 4", arg, sep=""))
   if(strands=="both") {
      system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
   } else {
      system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
   }
}

# fix a broken HTML markup
system("sed 's/<\\/body><\\/html>//' seqs.txt.html > seqs.html")
