# TOOL norm-illumina-SNP.R: "Illumina SNP arrays" (Illumina SNP array preprocessing. Input should be a tab-delimited text file with genotype calls. Typically such a file is created using GenCall software from Illumina.)
# INPUT chip.txt: chip.txt TYPE GENERIC 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER method: "Allel calling method" TYPE [illumina: illumina, crlmm: CRLMM] DEFAULT illumina (Allel calling method. Choosing Illumina, allels calls given in the input are accepted. For CRLMM, columns X Raw and Y Raw should be available in the file)
# PARAMETER cdfName: cdfName TYPE [default: default, human1mduov3b: human1mduov3b, human1mv1c: human1mv1c, human370quadv3c: human370quadv3c, human370v1c: human370v1c, human550v3b: human550v3b, human610quadv1b: human610quadv1b, human650v3a: human650v3a, human660quadv1a: human660quadv1a, humanomni1quadv1b: humanomni1quadv1b, humanomniexpress12v1b: humanomniexpress12v1b] DEFAULT default (Name of the description file to be used in normalization, by defualt obtained from data)
# PARAMETER SNRMin: "Signal-to-noise ratio" TYPE DECIMAL FROM 0 TO 10000 DEFAULT 5 (Value defining the minimum signal-to-noise ratio used to filter out samples, higher values mean more strict filtering. Affects only the CRLMM-function)

# Illumina SNP array normalization
# JTT 22.10.2008
# MK 21.05.2013 added crlmm normalization

# Reading data
if(method = "illumina") {
	#allel calls
	firstfield <- scan(dir(), what = "", sep = ",", flush = TRUE, quiet = TRUE, blank.lines.skip = FALSE, multi.line = FALSE)
	skip <- grep("[Data]", firstfield, fixed = TRUE)
	samples <- read.table(dir(), skip = skip, header = TRUE, sep = "\t", as.is = TRUE, check.names = FALSE, colClasses = "character")

	# Transforming the data to a wide data frame
	ids<- unique(samples$"Sample ID")
	dat2<-matrix(nrow=nrow(samples)/length(ids), ncol=length(ids), data=NA)
	for(i in 1:length(ids)) {
   		cursam<-samples[samples$"Sample ID"==ids[i],]
   		dat2[,i]<-as.numeric(as.factor(paste(cursam$"Allele1 - Top", cursam$"Allele2 - Top", sep="")))
	}
	dat2<-data.frame(dat2)
	rownames(dat2)<-cursam$"SNP Name"
	colnames(dat2)<-paste("chip.", ids, sep="")
} else {	
	library(crlmm)
	library(oligoClasses)
	library(human1mduov3bCrlmm)          #Illumina
	library(human1mv1cCrlmm)             #Illumina
	library(human370quadv3cCrlmm)        #Illumina
	library(human370v1cCrlmm)            #Illumina
	library(human550v3bCrlmm)            #Illumina
	library(human610quadv1bCrlmm)        #Illumina
	library(human650v3aCrlmm)            #Illumina
	library(human660quadv1aCrlmm)        #Illumina
	library(humanomni1quadv1bCrlmm)      #Illumina
	library(humanomniexpress12v1bCrlmm)  #Illumina
	
	if(cdfName == "default") {
		XY = readGenCallOutput(file="chip.txt", colnames=list("SampleID"="Sample ID", "SNPID"="SNP Name", "XRaw"="X Raw", "YRaw"="Y Raw"), type=list("SampleID"="character", "SNPID"="character", "XRaw"="integer", "YRaw"="integer"))
		crlmmOut <- try(crlmm::crlmmIllumina(XY=XY, verbose=FALSE, SNRMin=SNRMin))
	} else {
		XY = readGenCallOutput(file="chip.txt", cdfName=cdfName, colnames=list("SampleID"="Sample ID", "SNPID"="SNP Name", "XRaw"="X Raw", "YRaw"="Y Raw"), type=list("SampleID"="character", "SNPID"="character", "XRaw"="integer", "YRaw"="integer"))
		crlmmOut <- try(crlmm::crlmmIllumina(XY=XY, verbose=FALSE, cdfName=cdfName, SNRMin=SNRMin))
	}
	
	#Below commands allow to run crlmm on idat (scanner output) data 
	#samples <- read.csv("../idat/sampleSheet.csv", as.is=TRUE)			
	#readIdatFiles -method provide functionality to read idat files 
	#RG <- readIdatFiles(sampleSheet=samples, path= "../idat", arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"), fileExt=list(green="Grn.idat", red="Red.idat"))
	#crlmmOut <- try(crlmm::crlmmIllumina(RG=RG, verbose=FALSE, cdfName=cdfName, SNRMin=SNRMin))
				
	if(class(crlmmOut)=="try-error") {
		if(length(grep("No arrays above quality", crlmmOut[1])) > 0) {
			stop("CHIPSTER-NOTE: No arrays above quality threshold. Please try a lower signal-to-noise ratio parameter");
		} else {
			stop("CHIPSTER-NOTE: An unexpected error has occurred. Please contact Chipster support")
		}
	}
	
	genotypes<-calls(crlmmOut)
	flags<-callProbability(crlmmOut)
	flags2<-flags
	flags2[flags<=quantile(flags, 0.95)]<-"P"
	flags2[flags>quantile(flags, 0.95)]<-"M"
	colnames(genotypes)<-paste("chip.", colnames(genotypes), sep="")
	colnames(flags2)<-paste("flags.", colnames(flags), sep="")
	
	dat2 <- data.frame(genotypes, flags2)
}

# Writes out a phenodata table
chiptype<-"Illumina"
sample<-colnames(dat2)
group<-c(rep("", ncol(dat2)))
training<-c(rep("", ncol(dat2)))
time<-c(rep("", ncol(dat2)))
random<-c(rep("", ncol(dat2)))
write.table(data.frame(sample=sample, chiptype=chiptype, group=group, training=training), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Writing data to disk
write.table(dat2, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
