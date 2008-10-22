# ANALYSIS "Miscellaneous"/"Multiple sequence alignment" (Does a multiple sequence alignment to sequences. 
# This tool needs a sequence file in FastA format. Output is a sequence alignment in Clustal format.) 
# INPUT GENERIC seqs.txt OUTPUT aligned-seqs.txt 
# PARAMETER program [muscle] DEFAULT muscle (Which program to use)
# PARAMETER method [slow, fast] DEFAULT fast (Use slow or fast settings)
# PARAMETER seqtype [DNA, protein] DEFAULT DNA (Are the sequences DNA or protein)


# Sets the path the muscle executable
path.muscle<-c("/v/linux26_x86_64/appl/molbio/muscle/muscle3.6/muscle")

# Sets the path the mafft executable
path.mafft<-c("mafft")

# Runs the analysis
if(program=="muscle" & method=="slow" & seqtype=="DNA") {
   system(paste(path.muscle, " -in seqs.txt -out aligned-seqs.txt -clwstrict -seqtype dna -stable", sep=""))
}

if(program=="muscle" & method=="fast" & seqtype=="DNA") {
   system(paste(path.muscle, " -in seqs.txt -out aligned-seqs.txt -clwstrict -maxiters 1 -diags1 -seqtype dna -stable", sep=""))
}

if(program=="muscle" & method=="slow" & seqtype=="protein") {
   system(paste(path.muscle, " -in seqs.txt -out aligned-seqs.txt -clwstrict -seqtype protein -stable", sep=""))
}

if(program=="muscle" & method=="fast" & seqtype=="protein") {
   system(paste(path.muscle, " -in seqs.txt -out aligned-seqs.txt -clwstrict -seqtype protein -maxiters 1 -diags1 -sv -distance1 kbit20_3 -stable", sep=""))
}



if(program=="mafft" & method=="slow" & seqtype=="DNA") {
   write("module load mafft", "module.txt")
   write(paste(path.mafft, " --localpair --maxiterate 1000 --clustalout --inputorder --nuc seqs.txt > aligned-seqs.txt", sep=""), "module.txt", append=T)
   system("tcsh module.txt")
}

if(program=="mafft" & method=="slow" & seqtype=="protein") {
   write("module load mafft", "module.txt")
   write(paste(path.mafft, " --localpair --maxiterate 1000 --clustalout --inputorder --amino seqs.txt > aligned-seqs.txt", sep=""), "module.txt", append=T)
   system("tcsh module.txt")
}

if(program=="mafft" & method=="fast" & seqtype=="DNA") {
   write("module load mafft", "module.txt")
   write(paste(path.mafft, " --retree 2 --maxiterate 0 --clustalout --inputorder --nuc seqs.txt > aligned-seqs.txt", sep=""), "module.txt", append=T)
   system("tcsh module.txt")
}

if(program=="mafft" & method=="fast" & seqtype=="protein") {
   write("module load mafft", "module.txt")
   write(paste(path.mafft, " --retree 2 --maxiterate 0 --clustalout --inputorder --amino seqs.txt > aligned-seqs.txt", sep=""), "module.txt", append=T)
   system("tcsh module.txt")
}

if(program=="clustal" & method=="slow") {
   library(Biostrings)
   library(dna)
   fas<-readFASTA("seqs.txt")
   v<-rep(NA, length(fas))
   for(i in 1:length(fas)) {
      v[i]<-fas[[i]]$desc
   }
   s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
   for(i in 1:length(fas)) {
      s[i,]<-fas[[i]]$seq
   }
   names(s)<-"Sequence"
   v<-gsub(">", "", v)
   rownames(s)<-v
   alignment<-clustalw(as.matrix(s), names=rownames(s))
}

if(program=="clustal" & method=="fast") {
   library(Biostrings)
   library(dna)
   fas<-readFASTA("seqs.txt")
   v<-rep(NA, length(fas))
   for(i in 1:length(fas)) {
      v[i]<-fas[[i]]$desc
   }
   s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
   for(i in 1:length(fas)) {
      s[i,]<-fas[[i]]$seq
   }
   names(s)<-"Sequence"
   v<-gsub(">", "", v)
   rownames(s)<-v
   alignment<-clustalw(as.matrix(s), names=rownames(s), quick.pairalign=T)
}

if(program=="clustal") {
   write(file="aligned-seqs.txt", "CLUSTAL W (1.83) multiple sequence alignment")
   write(file="aligned-seqs.txt", cat("\n"), append=T)
   for(i in 1:nrow(alignment$seqout)) {
      if(i<nrow(alignment$seqout)) {
         write(file="aligned-seqs.txt", x=paste(rownames(alignment$seqout)[i], alignment$seqout[i,1], sep="\t"), append=T) 
      } else {
         write(file="aligned-seqs.txt", x=paste("", alignment$seqout[i,1], sep="\t"), append=T) 
      }    
   }
}
