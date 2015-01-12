# TOOL norm-affy-snp.R: "Affymetrix SNP arrays" (Affymetrix SNP array preprocessing using CEL-files. Probe sets are automatically flagged using P M flags.)
# INPUT microarray{...}.cel: microarray{...}.cel TYPE AFFY 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER version: CRLMM-version TYPE [oligo: oligo, crlmm: crlmm] DEFAULT crlmm (Use CRLMM function from the package oligo, i.e. the old version, or from the package crlmm, i.e. the new version)
# PARAMETER cdfName: "CDF name" TYPE [default: default, genomewidesnp5: genomewidesnp5, genomewidesnp6: genomewidesnp6] DEFAULT default (Name of the description file to be used in normalization, by default obtained from the CEL files. Affects only the function of the crlmm-package)
# PARAMETER SNRMin: "Signal-to-noise ratio" TYPE DECIMAL FROM 0 TO 10000 DEFAULT 5 (Value defining the minimum signal-to-noise ratio used to filter out samples, higher values mean more strict filtering. Affects only the function of the crlmm-package)

# JTT: 24.4.2008 created
# MG: 26.10.2009 to cope with Affymetrix Human SNP 5.0 and 6.0 arrays
# MK: 17.05.2013 added option to use the new

# setwd(paste(getwd(), "/affy-snp", sep=""))

# Loading the libraries
if(version == "oligo") {
	library(oligo)
	library(oligoClasses)
	library(affyio)
} else {
	library(crlmm)
	library(oligoClasses)
	library(affyio)
}

library(pd.mapping50k.hind240)
library(pd.mapping50k.xba240)
library(pd.mapping250k.nsp)
library(pd.mapping250k.sty)
library(pd.genomewidesnp.5)
library(pd.genomewidesnp.6)
library(genomewidesnp5Crlmm)
library(genomewidesnp6Crlmm)

# Setting up the path to the data files
fullFilenames <- list.celfiles(path = getwd(), full.names = TRUE)

# Calculating the genotype calls
if(version == "oligo") {
	this.cdfName <- gsub("_", "", tolower(read.celfile.header(fullFilenames[1])$cdfName))
	if(this.cdfName != "genomewidesnp5" && this.cdfName != "genomewidesnp6") {
		stop("CHIPSTER-NOTE: Chip type not supported by oligo.")
	}	

	crlmmOut <- try(crlmm(fullFilenames, file.path(getwd(), "crlmmResults")))
	crlmmOut <- try(getCrlmmSummaries(file.path(getwd(), "crlmmResults")))
} else {
	if(cdfName != "default") {
		crlmmOut <- try(crlmm::crlmm(fullFilenames, verbose=FALSE, cdfName=cdfName, SNRMin=SNRMin))
	} else {	
		this.cdfName <- gsub("_", "", tolower(read.celfile.header(fullFilenames[1])$cdfName))
		if(this.cdfName != "genomewidesnp5" && this.cdfName != "genomewidesnp6") {
			stop("CHIPSTER-NOTE: Chip type not supported by CRLMM. Please choose a suitable description file from the drop-down box or use another pre-processing method")
		}	
			
		crlmmOut <- try(crlmm::crlmm(fullFilenames, verbose=FALSE, SNRMin=SNRMin))
	}
}

if(class(crlmmOut)=="try-error") {
	if(length(grep("No arrays above quality", crlmmOut[1])) > 0) {
		stop("CHIPSTER-NOTE: No arrays above quality threshold. Please try a lower signal-to-noise ratio parameter");
	} else {
		stop("CHIPSTER-NOTE: An unexpected error has occurred. Please contact Chipster support")
	}
}
	
# Putting the data in a correct format
genotypes<-calls(crlmmOut)
flags<-confs(crlmmOut)
flags2<-flags
flags2[flags<quantile(flags, 0.95, na.rm=T)]<-"P"
flags2[flags>=quantile(flags, 0.95, na.rm=T)]<-"M"
colnames(genotypes)<-paste("chip.", colnames(genotypes), sep="")
colnames(flags2)<-paste("flags.", colnames(flags), sep="")

# Writes out a phenodata table
chiptype<-c(rep("affy-SNP", ncol(genotypes)))
sample<-colnames(genotypes)
group<-c(rep("", ncol(genotypes)))
write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Writing the results into files
write.table(data.frame(genotypes, flags2), "normalized.tsv", row.names=T, col.names=T, quote=F, sep="\t")
