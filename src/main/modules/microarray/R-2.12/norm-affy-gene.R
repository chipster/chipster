# TOOL norm-affy-gene.R: "Affymetrix gene arrays" (Affymetrix RMA preprocessing for CEL-files. Please note that the preprocessing might take a long time. YOU HAVE TO SPECIFY THE CHIPTYPE.)
# INPUT microarray{...}.cel: microarray{...}.cel TYPE AFFY 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER chiptype: chiptype TYPE [empty: empty, human-1.0-ST: human-1.0-ST, human-1.1-ST: human-1.1-ST, mouse-1.0-ST: mouse-1.0-ST, mouse-1.1-ST: mouse-1.1-ST, rat-1.0-ST: rat-1.0-ST, rat-1.1-ST: rat-1.1-ST] DEFAULT empty (Chiptype)


# Affymetrix normalization
# JTT, 3.2.2009
# MG. 18.10.2011, added version 1.1 ST arrays

# Initializes analyses
library(affy)

# Reads in data
dat<-ReadAffy()
if(chiptype=="empty") {
	stop("CHIPSTER-NOTE: You need to specify the chiptype. Please run the tool again.")
}

if(chiptype=="human-1.0-ST") {
	dat@cdfName<-"hugene10stv1hsentrezgcdf"
	dat@annotation<-"hugene10stv1hsentrezgcdf"
	chiptype<-"hugene10stv1hsentrezg.db"
}
if(chiptype=="human-1.1-ST") {
	dat@cdfName<-"hugene11stv1hsentrezgcdf"
	dat@annotation<-"hugene11stv1hsentrezgcdf"
	chiptype<-"hugene11stv1hsentrezg.db"
}
if(chiptype=="mouse-1.0-ST") {
	dat@cdfName<-"mogene10stv1mmentrezgcdf"
	dat@annotation<-"mogene10stv1mmentrezgcdf"
	chiptype<-"mogene10stv1mmentrezg.db"
}
if(chiptype=="mouse-1.1-ST") {
	dat@cdfName<-"mogene11stv1mmentrezgcdf"
	dat@annotation<-"mogene11stv1mmentrezgcdf"
	chiptype<-"mogene11stv1mmentrezg.db"
}
if(chiptype=="rat-1.0-ST") {
	dat@cdfName<-"ragene10stv1rnentrezgcdf"
	dat@annotation<-"ragene10stv1rnentrezgcdf"
	chiptype<-"ragene10stv1rnentrezg.db"
}
if(chiptype=="rat-1.1-ST") {
	dat@cdfName<-"ragene11stv1rnentrezgcdf"
	dat@annotation<-"ragene11stv1rnentrezgcdf"
	chiptype<-"ragene11stv1rnentrezg.db"
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

# Writing out data
a<-try(library(chiptype, character.only=T))
if(chiptype!="empty" & class(a)!="try-error") {
	# Including gene names to data
	lib2<-sub('.db','',chiptype)
	symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep="")))))[rownames(dat2),])
	genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep="")))))[rownames(dat2),])
	# Fxes an issue introduced in BioC2.4 where the "#" character is introduced in some gene names
	genename <- gsub("#", "", genename)
	symbol <- gsub("'", "", symbol)
	genename <- gsub("'", "", genename)
	# Writes the results into a file
	write.table(data.frame(symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
} 

if(chiptype=="empty" | class(a)=="try-error") {
	write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
}
