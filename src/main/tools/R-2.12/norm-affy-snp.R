# ANALYSIS Normalisation/"Affymetrix SNP arrays" (Affymetrix SNP array preprocessing using CEL-files. 
# Probe sets are automatically flagged using P/M flags.) 
# INPUT AFFY microarray[...].cel OUTPUT normalized.tsv, phenodata.tsv 


# Affymetrix SNP array normalization
# JTT 24.4.2008
#
# modified by MG, 26.10.2009
# to cope with Affymetrix Human SNP 5.0 and 6.0 arrays

# setwd(paste(getwd(), "/affy-snp", sep=""))

# Loading the libraries
library(oligo)
library(pd.mapping50k.hind240)
library(pd.mapping50k.xba240)
library(pd.mapping250k.nsp)
library(pd.mapping250k.sty)
library(pd.genomewidesnp.5)
library(pd.genomewidesnp.6)

# Setting up the path to the data files
fullFilenames <- list.celfiles(path = getwd(), full.names = TRUE)

# Calculating the genotype calls
crlmm(fullFilenames, file.path(getwd(), "crlmmResults"))
crlmmOut <- getCrlmmSummaries(file.path(getwd(), "crlmmResults"))

# Putting the data in a correct format
genotypes<-calls(crlmmOut)
flags<-confs(crlmmOut)
flags2<-flags
flags2[flags<=quantile(flags, 0.95)]<-"P"
flags2[flags>quantile(flags, 0.95)]<-"M"
colnames(genotypes)<-paste("chip.", colnames(genotypes), sep="")
colnames(flags2)<-paste("flags.", colnames(flags), sep="")

# Writes out a phenodata table
chiptype<-c(rep("affy-SNP", ncol(genotypes)))
sample<-colnames(genotypes)
group<-c(rep("", ncol(genotypes)))
write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Writing the results into files
write.table(data.frame(genotypes, flags2), "normalized.tsv", row.names=T, col.names=T, quote=F, sep="\t")
