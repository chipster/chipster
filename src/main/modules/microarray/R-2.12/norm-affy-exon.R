# TOOL norm-affy-exon.R: "Affymetrix exon arrays" (Affymetrix RMA preprocessing for CEL-files. YOU HAVE TO SPECIFY THE CHIPTYPE.)
# INPUT microarray{...}.cel: microarray{...}.cel TYPE AFFY 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER chiptype: chiptype TYPE [empty: empty, human: human, mouse: mouse, rat: rat] DEFAULT empty ()
# PARAMETER summary.feature: summary.feature TYPE [gene: gene, exon: exon] DEFAULT gene (Output summary type)


# Affymetrix normalization
# JTT 8.6.2006
# Changes to column naming on 29.6.2006
# Changes to phenodata table writing on 29.1.2007
#
# modified 11.11.2009, MG
# Changes to CDF package names to match BRainArray version 12
# 
# modified 25.11.2009, MG
# Gene symbols and gene names are now incorporated to the table
# of normalized values when using "gene" as summary feature
# but not for "exon", since the annotation packages do not
# support this

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
	dat@cdfName<-"huex10stv2hsentrezgcdf"
	dat@annotation<-"huex10stv2hsentrezg.db"
}
if(chiptype=="mouse" & summary.feature=="gene") {
	dat@cdfName<-"moex10stv1mmentrezgcdf"
	dat@annotation<-"moex10stv1mmentrezg.db"
}
if(chiptype=="rat" & summary.feature=="gene") {
	dat@cdfName<-"raex10stv1rnentrezgcdf"
	dat@annotation<-"raex10stv1rnentrezg.db"
}
chiptype<-dat@annotation

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

# Writing out data
a<-try(library(chiptype, character.only=T))
if(chiptype!="empty" & class(a)!="try-error") {
	if (summary.feature=="exon") {
		# Writes the results into a file
		write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	}
	if (summary.feature=="gene") {
		
		# Including gene names to data
		lib2<-sub('.db','',chiptype)
		symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep="")))))[rownames(dat2),])
		genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep="")))))[rownames(dat2),])
		
		# Fixes an issue introduced in BioC2.4 where the "#" character is introduced in some gene names
		genename <- gsub("#", "", genename)
		
		# Writes the results into a file
		write.table(data.frame(symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	}
} 

