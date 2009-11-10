# ANALYSIS Normalisation/Affymetrix (Affymetrix preprocessing using CEL-files. Probe sets are automatically flagged 
# using P/A/M flags. Variance stabilization can be applied only with MAS5 or Plier preprocessing methods. Custom
# chiptype can't be used with Plier preprocessing.)
# INPUT AFFY microarray[...].cel OUTPUT normalized.tsv, phenodata.tsv 
# PARAMETER normalization.method [mas5, plier, rma, gcrma, li-wong] DEFAULT rma (Preprocessing method)
# PARAMETER stabilize.variance [yes, no] DEFAULT no (Variance stabilazing normalization)
# PARAMETER custom.chiptype [empty, hs133ahsentrezg (hgu133a), hs133av2hsentrezg (hgu133av2), hs133phsentrezg (hgu133plus2), hs133bhsentrezg (hgu133b), hs95av2hsentrezg (hgu95av2), mm430ammentrezg (moe430), mm430a2mmentrezg (moe4302), mm430bmmentrezg (moe430b), mm430mmentrezg (moe430plus2), mm74av1mmentrezg (mgu74a), mm74av2mmentrezg (mgu74av2), mm74bv2mmentrezg (mgu74bv2), mm74cv2mmentrezg (mgu74cv2), rn230arnensg (rgu230a), rn230brnentrezg (rgu230b), rn230rnentrezg (rat2302), rn34arnentrezg (rgu34a)] DEFAULT empty (custom chiptype)


# Affymetrix normalization
# JTT 8.6.2006
# Changes to column naming on 29.6.2006
# Changes to phenodata table writing on 29.1.2007

# Renaming variables
norm<-normalization.method
stabvar<-stabilize.variance

# Initializes analyses
library(affy)
library(gcrma)

# Reads in data
dat<-ReadAffy()

# Modifies the data objects to take custom chiptype into account
if(custom.chiptype!="empty") {
   chiptype<-custom.chiptype
   library(Biostrings)
   chiptype<-substr(x=chiptype, start=1, stop=(matchPattern("(", chiptype)@start-1))
   dat@annotation<-chiptype
   dat@cdfName<-chiptype
} else {
   chiptype<-dat@annotation
}

# Normalizations (MAS5, RMA, GCRMA or Li-Wong (dChip))

# MAS5 normalization
if(norm=="mas5") {
   dat2<-exprs(mas5(dat))
   calls<-exprs(mas5calls(dat))
   if(stabvar=="yes") {
      if(ncol(dat2)<2) {
         stop("You need to have at least two chip to be able to use VSN!")
      }
      library(vsn)
      dat2<-exprs(vsn(dat2))
   } else {
      dat2<-log2(dat2)
   }
   dat2<-as.data.frame(round(dat2, digits=2))
   calls<-as.data.frame(calls)
   names(dat2)<-paste("chip.", names(dat2), sep="")
   names(calls)<-paste("flag.", names(calls), sep="")
   dat2<-data.frame(dat2, calls)
}

# PLIER normalization
if(norm=="plier" & custom.chiptype=="empty") {
   library(plier)
   dat2<-exprs(justPlier(eset=dat,replicate=1:length(dat),get.affinities=FALSE,normalize=FALSE,norm.type=c("together"),augmentation=0.1,defaultaffinity=1.0,defaultconcentration=1.0,attenuation=0.005,seaconvergence=0.000001,seaiteration=3000,gmcutoff=0.15,probepenalty=0.001,concpenalty=0.000001,usemm=TRUE,usemodel=FALSE,fitaffinity=T,plierconvergence=0.000001,plieriteration=3000,dropmax=3.0,lambdalimit=0.01,optimization=0))
   calls<-exprs(mas5calls(dat))
   if(stabvar=="yes") {
      if(ncol(dat2)<2) {
         stop("You need to have at least two chip to be able to use VSN!")
      }
      library(vsn)
      dat2<-exprs(vsn(dat2))
   } else {
      dat2<-log2(dat2)
   }
   dat2<-as.data.frame(round(dat2, digits=2))
   calls<-as.data.frame(calls)
   names(dat2)<-paste("chip.", names(dat2), sep="")
   names(calls)<-paste("flag.", names(calls), sep="")
   dat2<-data.frame(dat2, calls)
}
if(norm=="plier" & custom.chiptype!="empty") {
   stop("Custom chipstypes can't be used with Plier! Use some other preprocessing method.")
}

# RMA normalization
if(norm=="rma" & custom.chiptype=="empty") {
   dat2<-exprs(justRMA())
   calls<-exprs(mas5calls(dat))
   dat2<-as.data.frame(round(dat2, digits=2))
   calls<-as.data.frame(calls)
   names(dat2)<-paste("chip.", names(dat2), sep="")
   names(calls)<-paste("flag.", names(calls), sep="")
   dat2<-data.frame(dat2, calls)
}
if(norm=="rma" & custom.chiptype!="empty") {
   dat2<-exprs(rma(dat))
   calls<-exprs(mas5calls(dat))
   dat2<-as.data.frame(round(dat2, digits=2))
   calls<-as.data.frame(calls)
   names(dat2)<-paste("chip.", names(dat2), sep="")
   names(calls)<-paste("flag.", names(calls), sep="")
   dat2<-data.frame(dat2, calls)
}


# GCRMA normalization
if(norm=="gcrma" & custom.chiptype=="empty") {
   dat2<-exprs(justGCRMA(type=c("fullmodel"), fast=T, optimize.by=c("speed")))
   calls<-exprs(mas5calls(dat))
   dat2<-as.data.frame(round(dat2, digits=2))
   calls<-as.data.frame(calls)
   names(dat2)<-paste("chip.", names(dat2), sep="")
   names(calls)<-paste("flag.", names(calls), sep="")
   dat2<-data.frame(dat2, calls)
}
if(norm=="gcrma" & custom.chiptype!="empty") {
   dat2<-exprs(gcrma(dat))
   calls<-exprs(mas5calls(dat))
   dat2<-as.data.frame(round(dat2, digits=2))
   calls<-as.data.frame(calls)
   names(dat2)<-paste("chip.", names(dat2), sep="")
   names(calls)<-paste("flag.", names(calls), sep="")
   dat2<-data.frame(dat2, calls)
}


# Li-Wong (dChip) normalization
if(norm=="li-wong") {
   dat2<-exprs(expresso(dat, normalize.method="invariantset", bg.correct=FALSE, pmcorrect.method="pmonly", summary.method="liwong"))
   calls<-exprs(mas5calls(dat))
   dat2<-as.data.frame(round(dat2, digits=2))
   calls<-as.data.frame(calls)
   names(dat2)<-paste("chip.", names(dat2), sep="")
   names(calls)<-paste("flag.", names(calls), sep="")
   dat2<-data.frame(dat2, calls)
}

# Writes out a phenodata table
sample<-rownames(pData(dat))
group<-c(rep("", nrow(pData(dat))))
training<-c(rep("", nrow(pData(dat))))
time<-c(rep("", nrow(pData(dat))))
random<-c(rep("", nrow(pData(dat))))
write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Writing out data
a<-try(library(chiptype, character.only=T))
if(chiptype!="empty" & class(a)!="try-error") {
   # Including gene names to data
   symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(chiptype, "SYMBOL", sep="")))))[rownames(dat2),])
   genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(chiptype, "GENENAME", sep="")))))[rownames(dat2),])
   # Writes the results into a file
   write.table(data.frame(symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
} 

if(chiptype=="empty" | class(a)=="try-error") {
   write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
}
