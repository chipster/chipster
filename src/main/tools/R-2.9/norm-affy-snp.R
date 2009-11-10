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

# Calculating the intensity values for the probes
# intval<-justSNPRMA(filenames=fullFilenames, verbose=T, normalizeToHapmap=T)

# Calculating the genotype calls
aboutSamples <- data.frame(gender = c(rep("male", length(fullFilenames))))
rownames(aboutSamples) <- basename(fullFilenames)
aboutVars <- data.frame(labelDescription = "male/female")
rownames(aboutVars) <- "gender"
pd <- new("AnnotatedDataFrame", data = aboutSamples, varMetadata = aboutVars)
crl <- justCRLMM(fullFilenames, phenoData=pd, verbose=T)

# Putting the data in a correct format
crl2<-as.list(crl@assayData)
genotypes<-crl2$calls
flags<-crl2$callsConfidence
for(i in 1:nrow(flags)) {
   flags[i,]<-ifelse(flags[i,]<=0.95, "M", "P")
}
colnames(genotypes)<-paste("chip.", colnames(genotypes), sep="")
colnames(flags)<-paste("flags.", colnames(flags), sep="")

# Writes out a phenodata table
chiptype<-crl@annotation
sample<-rownames(pData(crl))
group<-c(rep("", nrow(pData(crl))))
training<-c(rep("", nrow(pData(crl))))
time<-c(rep("", nrow(pData(crl))))
random<-c(rep("", nrow(pData(crl))))
write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Writing the results into files
write.table(data.frame(genotypes, flags), "normalized.tsv", row.names=T, col.names=T, quote=F, sep="\t")
