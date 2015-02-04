# TOOL norm-affy-exon.R: "Affymetrix exon arrays" (Affymetrix RMA preprocessing for CEL-files. YOU HAVE TO SPECIFY THE CHIPTYPE.)
# INPUT microarray{...}.cel: microarray{...}.cel TYPE AFFY 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER chiptype: "Chiptype" TYPE [empty: empty, human: human, human-hta20: human-hta20, mouse: mouse, rat: rat] DEFAULT empty (Chiptype)
# PARAMETER summary.feature: "Summary feature" TYPE [gene: gene, exon: exon] DEFAULT gene (Output summary type)

# JTT 08.06.2006: Created
# JTT 29.06.2006: Changes to column naming on 
# JTT 29.01.2007: Changes to phenodata table writing on 
# MG 11.11.2009: modified 11.11.2009
# MG 25.11.2009: modified 25.11.2009
# MK 12.02.2014: added hta-20 chips to the body.

# Initializes analyses
library(oligo)
library(affy)
library(biomaRt)

# Reads in data
if(chiptype=="empty") {
	stop("You need to specify the chiptype. Please run the script again.")
}
if(chiptype=="human" & summary.feature=="exon") {
	custom_cdf = "huex10sthsensecdf"
	chiptype <- "huex10sthsensecdf.db"	
	#dat@cdfName<-"huex10stv2hsensecdf"
	#dat@annotation<-"huex10stv2hsensecdf.db"
}
if(chiptype=="mouse" & summary.feature=="exon") {
	custom_cdf <- "moex10stmmensecdf"
	chiptype <- "moex10stmmensecdf.db"	
	#dat@cdfName<-"moex10stv1mmensecdf"
	#dat@annotation<-"moex10stv1mmensecdf.db"
}
if(chiptype=="rat" & summary.feature=="exon") {
	custom_cdf <- "raex10strnensecdf"
	chiptype <- "raex10strnensecdf.db"	
	#dat@cdfName<-"raex10stv1rnensecdf"
	#dat@annotation<-"raex10stv1rnensecdf.db"
}
if(chiptype=="human-hta20" & summary.feature=="exon") {
	custom_cdf = "hta20hsensecdf"
	chiptype <- "hta20hsensecdf.db"	
	#dat@cdfName<-"huex10stv2hsensecdf"
	#dat@annotation<-"huex10stv2hsensecdf.db"
}

if(chiptype=="human" & summary.feature=="gene") {
	custom_cdf <- "huex10sthsentrezgcdf"
	chiptype <- "huex10sthsentrezg.db"	
	#dat@cdfName<-"huex10stv2hsentrezgcdf"
	#dat@annotation<-"huex10stv2hsentrezg.db"
}
if(chiptype=="mouse" & summary.feature=="gene") {
	custom_cdf <- "moex10stmmentrezgcdf"
	chiptype <- "moex10stmmentrezg.db"
	#dat@cdfName<-"moex10stv1mmentrezgcdf"
	#dat@annotation<-"moex10stv1mmentrezg.db"
}
if(chiptype=="rat" & summary.feature=="gene") {
	custom_cdf <- "raex10strnentrezgcdf"
	chiptype <- "raex10strnentrezg.db"
	#dat@cdfName<-"raex10stv1rnentrezgcdf"
	#dat@annotation<-"raex10stv1rnentrezg.db"
}
if(chiptype=="human-hta20" & summary.feature=="gene") {
	custom_cdf = "hta20hsentrezgcdf"
	chiptype <- "hta20hsentrezg.db"	
	#dat@cdfName<-"huex10stv2hsensecdf"
	#dat@annotation<-"huex10stv2hsensecdf.db"
}

#chiptype<-dat@annotation

# Normalizations
if(chiptype != "oligo") {
	dat2 <- justRMA(filenames=list.celfiles(), cdfname=custom_cdf)
} else {
	data.raw <- read.celfiles(filenames=list.celfiles())
	dat2 <- oligo::rma(data.raw)
}

dat2<-as.data.frame(round(exprs(dat2), digits=2))
names(dat2)<-paste("chip.", names(dat2), sep="")

# Writes out a phenodata
sample<-gsub("chip.", "", colnames(dat2))
group<-c(rep("", ncol(dat2)))
training<-c(rep("", ncol(dat2)))
time<-c(rep("", ncol(dat2)))
random<-c(rep("", ncol(dat2)))
chiptype<-chiptype
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
		symbol <- gsub("'", "", symbol)
		genename <- gsub("'", "", genename)
		
		# Writes the results into a file
		write.table(data.frame(symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	}
} else {
		write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)	
}

