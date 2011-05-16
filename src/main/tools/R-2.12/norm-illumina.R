# ANALYSIS Normalisation/"Illumina" (Illumina preprocessing using individual files. Every file includes data for one
# array, i.e., the data that has been imported through Import tool. YOU HAVE TO SPECIFY THE CHIPTYPE.)
# INPUT CDNA microarray[...].tsv OUTPUT normalized.tsv, phenodata.tsv
# PARAMETER normalize.chips [none, scale, quantile, vsn] DEFAULT quantile (Between arrays normalization method)
# PARAMETER beadstudio.version [1, 2, 3] DEFAULT 1 (BeadStudio version number) 
# PARAMETER chiptype [empty, Human-6v1, HumanRef-8v1, Human-6v2, HumanRef-8v2, Human-6v3, HumanRef-8v3, Human-HT12, Human-HT12v4, Mouse-6v1.0a, MouseRef-8v1.0a, Mouse-6v1.1, MouseRef-8v1.1, Mouse-6v2, MouseRef-8v2, RatRef-12] DEFAULT empty (chiptype)
# PARAMETER id.type [TargetID, ProbeID] DEFAULT TargetID (Which annotations to use) 
# PARAMETER produce.flags [yes, no] DEFAULT no (Automatic recoding of Detection-value as flags)

# Illumina data preprocessing and normalization for separate files
# JTT 17.10.2007

#normalize.chips<-"quantile"
#beadstudio.version<-c(1)
#chiptype<-"Human-6v1"
#id.type<-"TargetID"
#produce.flags<-"no"

# PARAMETER produce.flags [yes, no] DEFAULT no (Automatic recoding of Detection-value as flags)
# produce.flags<-c("no")

# Loads the libraries
library(limma)

# Renaming variables
normba<-normalize.chips

# Reading data
columns<-list(R="sample", Rb="sample", G="sample", Gb="sample")
annotation<-c("identifier")
columns.other<-c("flag")

files<-dir()
files<-files[files!="phenodata.tsv"]
dat<-read.maimages(files=files, columns=columns, annotation=annotation, other.columns=columns.other) 

# Normalization
dat2<-normalizeBetweenArrays(dat$R, method=normba)
if(normba!="vsn") {
   dat2<-log2(dat2)
} else {
   dat2<-dat2
}

# Rounding expression data to two digits
dat2<-round(dat2, digits=2)
dat2<-as.data.frame(dat2)
sample.names<-colnames(dat2)
sample.names<-paste(sample.names, ".tsv", sep="")
names(dat2)<-paste("chip.", sample.names, sep="")
rownames(dat2)<-dat$genes[,1]

# Producing flags
if(produce.flags=="yes") {
   if(length(dat$other$flag)!=0) {
      flags<-as.data.frame(dat$other$flag)
      flags2<-flags
      for(i in 1:ncol(flags)) {
         flags2[,i]<-as.numeric(as.vector(flags[,i]))
      }
      flags<-flags2
      names(flags)<-paste("flag.", names(flags), sep="")
   }
   if(length(dat$other$flag)==0) {
      flags<-matrix(nrow=0, ncol=0)
   }
}

if(produce.flags=="yes" & beadstudio.version==1) {
   m<-0.95
   a<-0.90
   flags[flags>m]<-"P"
   flags[flags>a & flags<=m]<-"M"
   flags[flags<=a]<-"A"
}

if(produce.flags=="yes" & beadstudio.version>1) {
   m<-0.05
   a<-0.10
   flags[flags>m]<-"P"
   flags[flags>a & flags<=m]<-"M"
   flags[flags<=a]<-"A"
}

# Writes out a phenodata table
group<-c(rep("", ncol(dat2)))
training<-c(rep("", ncol(dat2)))
time<-c(rep("", ncol(dat2)))
random<-c(rep("", ncol(dat2)))
if(id.type=="TargetID") {
   if(chiptype=="empty") {
      chiptype<-c("Illumina")
   }
   if(chiptype=="Human-6v1" | chiptype=="HumanRef-8v1") {
      chiptype<-c("illuminaHumanv1")
   }
   if(chiptype=="Human-6v2" | chiptype=="HumanRef-8v2") {
      chiptype<-c("illuminaHumanv2")
   }
   if(chiptype=="Human-6v3" | chiptype=="HumanRef-8v3") {
      chiptype<-c("illuminaHumanv3")
   }
   if(chiptype=="Human-HT12") {
      chiptype<-c("illuminaHumanv3")
   }
   if(chiptype=="Human-HT12v4") {
      chiptype<-c("illuminaHumanv4")
   }
   if(chiptype=="Mouse-6v1.0a" | chiptype=="MouseRef-8v1.0a") {
      chiptype<-c("illuminaMousev1")
   }
   if(chiptype=="Mouse-6v1.1" | chiptype=="MouseRef-8v1.1") {
      chiptype<-c("illuminaMousev1p1")
   }
   if(chiptype=="Mouse-6v2" | chiptype=="MouseRef-8v2") {
      chiptype<-c("illuminaMousev2")
   }
   if(chiptype=="RatRef-12") {
      chiptype<-c("illuminaRatv1")
   }
}

if(id.type=="ProbeID") {
   if(chiptype=="empty") {
      chiptype<-c("Illumina")
   }
   if(chiptype=="Human-6v1" | chiptype=="HumanRef-8v1") {
      chiptype<-c("illuminaHumanv1BeadID")
   }
   if(chiptype=="Human-6v2" | chiptype=="HumanRef-8v2") {
      chiptype<-c("illuminaHumanv2BeadID")
   }
   if(chiptype=="Human-6v3" | chiptype=="HumanRef-8v3") {
      chiptype<-c("illuminaHumanv3BeadID")
   }
   if(chiptype=="Human-HT12") {
      chiptype<-c("illuminaHumanv3BeadID")
   }
   if(chiptype=="Human-HT12v4") {
      chiptype<-c("illuminaHumanv4BeadID")
   }
   if(chiptype=="Mouse-6v1.0a" | chiptype=="MouseRef-8v1.0a") {
      chiptype<-c("illuminaMousev1BeadID")
   }
   if(chiptype=="Mouse-6v1.1" | chiptype=="MouseRef-8v1.1") {
      chiptype<-c("illuminaMousev1p1BeadID")
   }
   if(chiptype=="Mouse-6v2" | chiptype=="MouseRef-8v2") {
      chiptype<-c("illuminaMousev2BeadID")
   }
   if(chiptype=="RatRef-12") {
      chiptype<-c("illuminaRatv1BeadID")
   }
}
lib<-paste(chiptype, ".db", sep="")

# Write out a phenodata
write.table(data.frame(sample=sample.names, chiptype=lib, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

if(chiptype!="Illumina") {
   # Including gene names to data
   library(lib, character.only=T)
   symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(chiptype, "SYMBOL", sep="")))))[rownames(dat2),])
   genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(chiptype, "GENENAME", sep="")))))[rownames(dat2),])
   symbol<-gsub("#", "", symbol)
   genename<-gsub("#", "", genename)
   symbol <- gsub("'", "", symbol)
   genename <- gsub("'", "", genename)
   # Write out expression data
   if(produce.flags=="yes") {
      write.table(data.frame(symbol, description=genename, dat2, flags), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
   } else {
      write.table(data.frame(symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
   }
} else {
   if(produce.flags=="yes") {
      write.table(data.frame(dat2, flags), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
   } else {
      write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
   }
}
