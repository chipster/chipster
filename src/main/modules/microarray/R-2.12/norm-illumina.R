# TOOL norm-illumina.R: Illumina (Normalization of Illumina data. The data needs to be imported to Chipster using the Import tool, which producdes one file for each sample. YOU HAVE TO SPECIFY THE CHIPTYPE.)
# INPUT microarray{...}.tsv: microarray{...}.tsv TYPE CDNA 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER normalize.chips: "Normalization method" TYPE [none: none, scale: scale, quantile: quantile, vsn: vsn] DEFAULT quantile (Between arrays normalization method)
# PARAMETER beadstudio.version: "Illumina software version" TYPE [3: "GenomeStudio or BeadStudio 3", 2: "BeadStudio 2", 1: "BeadStudio 1"] DEFAULT 3 (Illumina software version)
# PARAMETER chiptype: "Chip type" TYPE [empty: empty, Human-6v1: Human-6v1, HumanRef-8v1: HumanRef-8v1, Human-6v2: Human-6v2, HumanRef-8v2: HumanRef-8v2, Human-6v3: Human-6v3, HumanRef-8v3: HumanRef-8v3, Human-HT12: Human-HT12, Human-HT12v4: Human-HT12v4, Mouse-6v1.0a: Mouse-6v1.0a, MouseRef-8v1.0a: MouseRef-8v1.0a, Mouse-6v1.1: Mouse-6v1.1, MouseRef-8v1.1: MouseRef-8v1.1, Mouse-6v2: Mouse-6v2, MouseRef-8v2: MouseRef-8v2, RatRef-12: RatRef-12] DEFAULT empty ()
# PARAMETER id.type: "Identifier type" TYPE [TargetID: TargetID, ProbeID: ProbeID] DEFAULT ProbeID (Which identifiers to use)
# PARAMETER OPTIONAL produce.flags: "Produce flags" TYPE [yes: yes, no: no] DEFAULT no (Automatic recording of Detection-value as flags)

# Illumina data preprocessing and normalization for separate files
# JTT 17.10.2007
# EK 3.4.2013 Changed parameter naming
# MK 20.06.2013 Stops if trying to produce flags from data that does not support this feature
# MK 26.05.2013 fixed bug in illumina detection p-value thresholds
# MK 27.05.2013 Illumina detection p-value cells marked as Ms are converted to 0.0 -values

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
      flags<-as.data.frame(dat$other$flag, stringsAsFactors=FALSE)
	  flags[flags=="M"] <- 0.0;
	  
	  flags2<-flags
      for(i in 1:ncol(flags)) {
         flags2[,i]<-as.numeric(as.vector(flags[,i]))
      }
      flags<-flags2
      names(flags)<-paste("flag.", names(flags), sep="")
  }
   if(length(dat$other$flag)==0) {
	  stop("CHIPSTER-NOTE: To produce detection values, your data need to contain flag columns representing detection p-value information")
      flags<-matrix(nrow=0, ncol=0)
   }
}

if(produce.flags=="yes" & beadstudio.version==1) {
	m<-0.95
   	a<-0.90
	flags_temp <- flags
   	flags[flags_temp>m]<-"P"
   	flags[flags_temp>a & flags_temp<=m]<-"M"
   	flags[flags_temp<=a]<-"A"
}

if(produce.flags=="yes" & beadstudio.version>1) {	
	m<-0.05
   	a<-0.10
	flags_temp <- flags;
	
	flags[flags_temp<m]<-"P"
	flags[flags_temp>=m & flags_temp<a]<-"M"
	flags[flags_temp>=a]<-"A"

	#flags[flags>m]<-"P"
   	#flags[flags>a & flags<=m]<-"M"
   	#flags[flags<=a]<-"A"
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
	   if(nrow(flags) == nrow(dat2)) {
           write.table(data.frame(symbol, description=genename, dat2, flags), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	   } else {
		   write.table(data.frame(symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	   }
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
