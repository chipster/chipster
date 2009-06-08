# ANALYSIS Normalisation/"Affymetrix gene arrays" (Affymetrix RMA preprocessing for CEL-files. 
# Please note that the preprocessing might take a long time, and a maximum of three human exon 
# arrays can be processed together. YOU HAVE TO SPECIFY THE CHIPTYPE.)
# INPUT AFFY microarray[...].cel OUTPUT normalized.tsv, phenodata.tsv 
# PARAMETER chiptype [empty, human, mouse, rat] DEFAULT empty (Chiptype)


# Affymetrix normalization
# JTT 3.2.2009

# Initializes analyses
library(affy)

# Reads in data
dat<-ReadAffy()
if(chiptype=="empty") {
   stop("You need to specify the chiptype. Please run the tool again.")
}

if(chiptype=="human") {
   dat@cdfName<-"hugene10stv1hsentrezgcdf"
   dat@annotation<-"hugene10stv1hsentrezgcdf"
   chiptype<-"hugene10stv1hsentrezg.db"
}
if(chiptype=="mouse") {
   dat@cdfName<-"mogene10stv1mmentrezgcdf"
   dat@annotation<-"mogene10stv1mmentrezgcdf"
   chiptype<-"mogene10stv1mmentrezg.db"
}
if(chiptype=="rat") {
   dat@cdfName<-"ragene10stv1rnentrezgcdf"
   dat@annotation<-"ragene10stv1rnentrezgcdf"
   chiptype<-"ragene10stv1rnentrezg.db"
}

# Normalizations
dat2<-exprs(rma(dat))
dat2<-as.data.frame(round(dat2, digits=2))
names(dat2)<-paste("chip.", names(dat2), sep="")

# Writes out a phenodata
sample<-rownames(pData(dat))
group<-c(rep("", nrow(pData(dat))))
training<-c(rep("", nrow(pData(dat))))
time<-c(rep("", nrow(pData(dat))))
random<-c(rep("", nrow(pData(dat))))
chiptype<-chiptype
write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Writes the results into a file
write.table(dat2, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
