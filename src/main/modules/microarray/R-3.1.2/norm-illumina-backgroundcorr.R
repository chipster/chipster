# TOOL norm-illumina-backgroundcorr.R: "Illumina with normexp background correction" (Normalization of Illumina data with normexp -background correction, using the nec-function. Two data files are needed: Summary probe profile file \(GenomeStudio output for each probe in the array\) and Summary control probe profile file \(GenomeStudio output file containing the control probe intensities\). YOU HAVE TO SPECIFY THE CHIPTYPE.)
# INPUT Sample_probe.txt: "Summary probe profile file" TYPE GENERIC
# INPUT Control_probe.txt: "Summary control probe profile file" TYPE GENERIC
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER normalize.chips: "Normalization method" TYPE [none: none, scale: scale, quantile: quantile] DEFAULT quantile (Between arrays normalization method)
# PARAMETER chiptype: "Chip type" TYPE [empty: empty, Human-6v1: Human-6v1, HumanRef-8v1: HumanRef-8v1, Human-6v2: Human-6v2, HumanRef-8v2: HumanRef-8v2, Human-6v3: Human-6v3, HumanRef-8v3: HumanRef-8v3, Human-HT12: Human-HT12, Human-HT12v4: Human-HT12v4, Mouse-6v1.0a: Mouse-6v1.0a, MouseRef-8v1.0a: MouseRef-8v1.0a, Mouse-6v1.1: Mouse-6v1.1, MouseRef-8v1.1: MouseRef-8v1.1, Mouse-6v2: Mouse-6v2, MouseRef-8v2: MouseRef-8v2, RatRef-12: RatRef-12] DEFAULT empty ()
# PARAMETER id.type: "Identifier type" TYPE [ProbeID: ProbeID] DEFAULT ProbeID (Which identifiers to use)
# PARAMETER OPTIONAL produce.flags: "Produce flags" TYPE [yes: yes, no: no] DEFAULT no (Automatic recording of Detection-value as flags)
# PARAMETER OPTIONAL annotations: "Include original annotations" TYPE [yes: yes, no: no] DEFAULT no (Include the original probe annotations. Note that these might be very outdated.)


# 10.4.2015 ML
# ei oo nyt eri beadstudioversioita.
# probeID:t ok, targetID -paketit puuttuu....
# vsn -paketti puuttui! vsn-osio poistettu.

library(limma)

dataFile   <- "Sample_probe.txt";  
qualFile   <- "Control_probe.txt";

data.raw <- read.ilmn(files=dataFile, ctrlfiles=qualFile, probeid="ProbeID")
# data.norm <- neqc(data.raw, negctrl="NEGATIVE") # neqc performs normexp background correction and quantile normalization aided by control probes

# Normalization (nec + normalizeBetweenArrays)
normba <- normalize.chips
data.norm <- nec(data.raw, negctrl="NEGATIVE")
if(normba!="vsn") {
	dat2 <- normalizeBetweenArrays(data.norm, method=normba)
} else {
	dat2 <- normalizeVSN(data.norm)
}
data.norm <- dat2

# Get the annotations. 
if(annotations=="yes") {
	orig_annotations <- data.norm$genes[,2]
}

# Detection values
if(produce.flags=="yes") {
	flags <- as.data.frame(data.norm$other)
}


# Phenodata -file

# rounding to 2 digits
dat2 <- data.norm$E
dat2<-as.data.frame(round(dat2, digits=2))

sample.names<-colnames(dat2)
#sample.names<-paste(sample.names, ".tsv", sep="")
names(dat2)<-paste("chip.", sample.names, sep="")


group<-c(rep("", ncol(dat2)))
training<-c(rep("", ncol(dat2)))
time<-c(rep("", ncol(dat2)))
random<-c(rep("", ncol(dat2)))

## which file type:
#if(id.type=="TargetID") {
#	if(chiptype=="empty") {
#		chiptype<-c("Illumina")
#	}
#	if(chiptype=="Human-6v1" | chiptype=="HumanRef-8v1") {
#		chiptype<-c("illuminaHumanv1")
#	}
#	if(chiptype=="Human-6v2" | chiptype=="HumanRef-8v2") {
#		chiptype<-c("illuminaHumanv2")
#	}
#	if(chiptype=="Human-6v3" | chiptype=="HumanRef-8v3") {
#		chiptype<-c("illuminaHumanv3")
#	}
#	if(chiptype=="Human-HT12") {
#		chiptype<-c("illuminaHumanv3")
#	}
#	if(chiptype=="Human-HT12v4") {
#		chiptype<-c("illuminaHumanv4")
#	}
#	if(chiptype=="Mouse-6v1.0a" | chiptype=="MouseRef-8v1.0a") {
#		chiptype<-c("illuminaMousev1")
#	}
#	if(chiptype=="Mouse-6v1.1" | chiptype=="MouseRef-8v1.1") {
#		chiptype<-c("illuminaMousev1p1")
#	}
#	if(chiptype=="Mouse-6v2" | chiptype=="MouseRef-8v2") {
#		chiptype<-c("illuminaMousev2")
#	}
#	if(chiptype=="RatRef-12") {
#		chiptype<-c("illuminaRatv1")
#	}
#}

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
	
	if(produce.flags=="yes" && annotations=="no") {
		write.table(data.frame(symbol, description=genename, dat2, flags), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=rownames(dat2))
	} else if(produce.flags=="yes" && annotations=="yes") {	   
		write.table(data.frame(symbol, description=genename, orig_annotations, dat2, flags), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=rownames(dat2))
	} else if(produce.flags=="no" && annotations=="yes") {	   
		write.table(data.frame(symbol, description=genename, orig_annotations, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=rownames(dat2))
	} else {
		write.table(data.frame(symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=rownames(dat2))
		
		
	}

} else {
	if(produce.flags=="yes" && annotations=="no") {
		write.table(data.frame( dat2, flags), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=rownames(dat2))
	} else if(produce.flags=="yes" && annotations=="yes") {  
		write.table(data.frame(orig_annotations, dat2, flags), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=rownames(dat2))
	} else if(produce.flags=="no" && annotations=="yes") {  
		write.table(data.frame(orig_annotations, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=rownames(dat2)) 
	} else {
		write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=rownames(dat2))
	}
}




