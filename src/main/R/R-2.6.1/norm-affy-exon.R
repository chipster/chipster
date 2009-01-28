# ANALYSIS Normalisation/"Affymetrix exon arrays" (Affymetrix RMA preprocessing for CEL-files. 
# Please note that the preprocessing might take a long time, and a maximum of three human exon 
# arrays can be processed together. YOU HAVE TO SPECIFY THE CHIPTYPE.)
# INPUT AFFY microarray[...].cel OUTPUT normalized.tsv, phenodata.tsv 
# PARAMETER chiptype [empty, human, mouse, rat] DEFAULT empty (Chiptype)
# PARAMETER summary.feature [gene, exon] DEFAULT gene (Output summary type)


# Affymetrix normalization
# JTT 8.6.2006
# Changes to column naming on 29.6.2006
# Changes to phenodata table writing on 29.1.2007

# Initializes analyses
library(affy)

# Reads in data
dat<-ReadAffy()
if(chiptype=="empty") {
   stop("You need to specify the chiptype. Please run the script again.")
}
if(chiptype=="human" & summary.feature=="exon") {
   dat@cdfName<-"exon.pmcdf"
   dat@annotation<-"exon.pmcdf"
}
if(chiptype=="mouse" & summary.feature=="exon") {
   dat@cdfName<-"mouseexonpmcdf"
   dat@annotation<-"mouseexonpmcdf"
}
if(chiptype=="rat" & summary.feature=="exon") {
   dat@cdfName<-"ratexonpmcdf"
   dat@annotation<-"ratexonpmcdf"
}

if(chiptype=="human" & summary.feature=="gene") {
   dat@cdfName<-"hsex10stv2hsentrezgcdf"
   dat@annotation<-"hsex10stv2hsentrezgcdf"
}
if(chiptype=="mouse" & summary.feature=="gene") {
   dat@cdfName<-"mmex10stv1mmentrezgcdf"
   dat@annotation<-"mmex10stv1mmentrezgcdf"
}
if(chiptype=="rat" & summary.feature=="gene") {
   dat@cdfName<-"rnex10stv1rnentrezgcdf"
   dat@annotation<-"rnex10stv1rnentrezgcdf"
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
chiptype<-dat@annotation
write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Writes the results into a file
write.table(dat2, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
