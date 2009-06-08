# ANALYSIS "Promoter Analysis"/"Known TFBSs" (Does a search for known transcription factor binding sites.)

date()

# Loads the libraries
library(Biostrings)

date()

# Reads in the FastA-formatted file containing the promoters
dat<-readFASTA("pro.txt", strip.descs=T, checkComments=T)

# Putting the sequences into a table
dat2<-matrix(ncol=2, nrow=length(dat), data=NA)
for(i in 1:length(dat)) {
   dat2[i,1]<- dat[[i]]$desc
   dat2[i,2]<- dat[[i]]$seq
}

# Reading in JASPARCORE matrices
mat<-read.table("CNE.txt", sep="\t", header=F)
matunique<-unique(mat[,1])

# Testing all JASPARCORE matrices against the promoter sequences
matfst<-matrix(nrow=length(matunique), ncol=length(dat), data=NA)
matrst<-matrix(nrow=length(matunique), ncol=length(dat), data=NA)
matfw<-matrix(nrow=length(matunique), ncol=length(dat), data=NA)
matrw<-matrix(nrow=length(matunique), ncol=length(dat), data=NA)
matfs<-matrix(nrow=length(matunique), ncol=length(dat), data=NA)
matrs<-matrix(nrow=length(matunique), ncol=length(dat), data=NA)

for(i in 1:length(matunique)) {
   mat2<-mat[mat[,1]==matunique[i],]
   pwm<-rbind(A=mat2[mat2[,2]=="A",][,4], C=mat2[mat2[,2]=="C",][,4], G=mat2[mat2[,2]=="G",][,4], T=mat2[mat2[,2]=="T",][,4])
   for(j in 1:nrow(dat2)) {
      seq<-DNAString(dat2[j,2])
      hitf<-matchPWM(pwm, seq)
      hitr<-matchPWM(reverseComplement(pwm), seq)

      hitfs<-rep(NA, length(hitf))
      hitrs<-rep(NA, length(hitr))
      if(length(hitf)>0) {
         for(k in 1:length(hitf)) {
            hitfs[k]<-as.character(hitf[[k]])
         }
      } 
      if(length(hitr)>0) {
         for(k in 1:length(hitr)) {
            hitrs[k]<-as.character(hitr[[k]])
         }
      }
      hitfst<-hitf@start
      hitrst<-hitr@start
      hitfw<-hitf@width
      hitrw<-hitf@width

      for(i in 1:length(hitfst)) {
         if(i==1) {
            hitfst2<-hitfst[i]
         } else {
            hitfst2<-paste(hitfst2, hitfst[i], sep=",")
         }   
      }

      for(i in 1:length(hitrst)) {
         if(i==1) {
            hitrst2<-hitrst[i]
         } else {
            hitrst2<-paste(hitrst2, hitrst[i], sep=",")
         }   
      }

      for(i in 1:length(hitfs)) {
         if(i==1) {
            hitfs2<-hitfs[i]
         } else {
            hitfs2<-paste(hitfs2, hitfs[i], sep=",")
         }   
      }

      for(i in 1:length(hitrs)) {
         if(i==1) {
            hitrs2<-hitrs[i]
         } else {
            hitrs2<-paste(hitrs2, hitrs[i], sep=",")
         }   
      }

      for(i in 1:length(hitfw)) {
         if(i==1) {
            hitfw2<-hitfw[i]
         } else {
            hitfw2<-paste(hitfw2, hitfw[i], sep=",")
         }   
      }

      for(i in 1:length(hitrw)) {
         if(i==1) {
            hitrw2<-hitrw[i]
         } else {
            hitrw2<-paste(hitrw2, hitrw[i], sep=",")
         }   
      }

      matfst[i,j]<-hitfst2
      matrst[i,j]<-hitrst2      
      matfw[i,j] <-hitfw2
      matrw[i,j] <-hitfw2 
      matfs[i,j] <-hitfs2
      matrs[i,j] <-hitrs2 

   }
}

# Now the matrices hold information for
# matfst & matrst		start positions
# matfs  & matfs        target sequences
# matfw  & matrw        target motif widths

# Visualizations

pdf(file="promoter.pdf", onefile=T, paper="a4r", width=11.7, height=8.3)

# Creating the main elements of the plot
plot(1, 1, type="n", xlim=c(1,1000), ylim=c(0,1), yaxt="n", bty="n", ylab="", xlab="")
yval<-0.8/ncol(matfst)
# Where to draw the lines
yvals<-yval*1:ncol(matfst)+0.1

# Adding the binding sites
# How many colors to generate?
nres<-rowSums(!is.na(matfst))
cols<-rainbow(sum(nres>0))

# Taking only the "significant" rows
matfst2<-matfst[nres>0,]
matrst2<-matrst[nres>0,]
matfs2<-matfs[nres>0,]
matrs2<-matrs[nres>0,]
matfw2<-matfw[nres>0,]
matrw2<-matrw[nres>0,]

# Replacing missing values with zeroes
matfst2[is.na(matfst2)]<-0
matrst2[is.na(matrst2)]<-0
matfs2[is.na(matfs2)]<-0
matrs2[is.na(matrs2)]<-0
matfw2[is.na(matfw2)]<-0
matrw2[is.na(matrw2)]<-0

# Matrices
mats<-as.vector(matunique[nres>0])

# Plotting the TFBS on the sequences
for(i in 1:ncol(matfst2)) {     # sekvenssien yli
   for(j in 1:nrow(matfst2)) {  # matriisien yli
      start<-as.numeric(strsplit(matfst2[j,i], ",")[[1]])/1 # Jos promoottorien pituus oli 1000bp, muutoin vaihda
      width<-as.numeric(strsplit(matfw2[j,i], ",")[[1]])/1   # Jos promoottorien pituus oli 1000bp, muutoin vaihda
      for(k in 1:length(start)) {
         polygon(x=c(start[k], start[k], start[k]+width[k], start[k]+width[k]), y=c(yvals[i], yvals[i]+yval/4, yvals[i]+yval/4, yvals[i]), col=cols[j], border=NA)
      }
      start<-as.numeric(strsplit(matrst2[j,i], ",")[[1]])/1  # Jos promoottorien pituus oli 1000bp, muutoin vaihda
      width<-as.numeric(strsplit(matrw2[j,i], ",")[[1]])/1  # Jos promoottorien pituus oli 1000bp, muutoin vaihda
      for(k in 1:length(start)) {
         polygon(x=c(start[k], start[k], start[k]+width[k], start[k]+width[k]), y=c(yvals[i], yvals[i]-yval/4, yvals[i]-yval/4, yvals[i]), col=cols[j], border=NA)
      }
   }
}

# Plotting the lines representing sequences
for(i in 1:ncol(matfst)) {
   lines(x=c(1,1000), y=c(yvals[i], yvals[i]))  
   mtext(text=dat2[i,1], side=2, at=yvals[i], las=1, cex=0.75) 
}

# Plotting the legend
for(i in 1:ncol(matfst2)) {
   text(x=-8, y=yvals[i]+yval/6, cex=0.75, labels="+")
   text(x=-8, y=yvals[i]-yval/6, cex=0.75, labels="-")
}
legend(x="top", fill=cols, legend=matunique[nres>0], ncol=min(6, length(cols)), cex=0.75, bty="n")
title(main="CNE-matrices")

dev.off()

date()


