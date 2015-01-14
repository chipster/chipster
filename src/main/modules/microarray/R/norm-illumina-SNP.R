# TOOL norm-illumina-SNP.R: "Illumina SNP arrays" (Illumina SNP array preprocessing. Input should be a tab-delimited text file with genotype calls. Typically such a file is created using GenCall software from Illumina or downloaded from a GEO.)
# INPUT chip.txt: chip.txt TYPE GENERIC 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER method: "Allel calling method" TYPE [illumina: illumina, crlmm: CRLMM] DEFAULT illumina (Allel calling method. Choosing Illumina, allels calls given in the input are accepted. For CRLMM, allell calls are computed from X-Raw and Y-Raw values)
# PARAMETER input: "Input format" TYPE [illumina: illumina, xy: xy-matrix] DEFAULT illumina (Input format. Illumina method requires output file created using GenCall software from Illumina where samples are listed one below another, whereas CRLMM also accepts xy-matrix.)
# PARAMETER cdfName: "CDF name" TYPE [default: default, human1mduov3b: human1mduov3b, human1mv1c: human1mv1c, human370quadv3c: human370quadv3c, human370v1c: human370v1c, human550v3b: human550v3b, human610quadv1b: human610quadv1b, human650v3a: human650v3a, human660quadv1a: human660quadv1a, humanomni1quadv1b: humanomni1quadv1b, humanomniexpress12v1b: humanomniexpress12v1b, humancytosnp12v2p1h: humancytosnp12v2p1h, humanomni25quadv1b: humanomni25quadv1b, humanomni5quadv1b: humanomni5quadv1b] DEFAULT default (Name of the description file to be used in normalization, by defualt obtained from data)
# PARAMETER SNRMin: "Signal-to-noise ratio" TYPE DECIMAL FROM 0 TO 10000 DEFAULT 5 (Value defining the minimum signal-to-noise ratio used to filter out samples, higher values mean more strict filtering. Affects only the CRLMM-function)

# JTT 22.10.2008: Illumina SNP array normalization
# MK  21.05.2013: added crlmm normalization. Add xy-matrix support

# Reading data
if(method == "illumina") {
 	if(input == "xy") {
 		stop("CHIPSTER-NOTE: Illumina type of preprocessing cannot be applied to xy-matrices");
 	}

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

	chiptype<-"Illumina-SNP"
	sample<-colnames(dat2)
	group<-c(rep("", ncol(dat2)))
	write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)
} else {	
	library(crlmm)
	library(Biobase)

	if(cdfName == "default") {
		stop("CHIPSTER-NOTE: please specify your chiptype");
	}

	if(input == "xy") {
		xy_mat <- read.table(file="chip.txt", sep="\t", header=T, comment.char="")
		snps <- xy_mat[,grep("Name", colnames(xy_mat), perl=T)]

		X <- xy_mat[, grep(".X.Raw", colnames(xy_mat), perl=T)]
		colnames(X) <- gsub(".X.Raw", "", colnames(X), perl=T)

		Y <- xy_mat[, grep(".Y.Raw", colnames(xy_mat), perl=T)]
		colnames(Y) <- gsub(".Y.Raw", "", colnames(Y), perl=T)
		Y <- Y[,match(colnames(X), colnames(Y))]

		if(min(table( c(colnames(X), colnames(Y)))) != 2) {
			stop("CHIPSTER-NOTE: Different number of raw X and Y columns found in the record");
		}

		zeroes <- matrix(0, ncol=ncol(X), nrow=nrow(X))
		zeroes = (X == "0" | Y == "0")
    	colnames(zeroes) <- colnames(Y) <- colnames(X)
    	rownames(X) <- rownames(Y) <- snps

		XY = new("NChannelSet", X = initializeBigMatrix(name = "X", nr = nrow(X), nc = ncol(X), vmode = "integer", initdata = as.matrix(X)),
								 Y = initializeBigMatrix(name = "Y", nr = nrow(Y), nc = ncol(Y), vmode = "integer", initdata = as.matrix(Y)),
								 zero = initializeBigMatrix(name = "zero", nr = nrow(zeroes), nc = ncol(zeroes), vmode = "integer", initdata = as.matrix(zeroes)), 
								 annotation = cdfName, storage.mode = "environment");
	} else {
		#overwrite the original function
		source(file.path(chipster.common.path, "crlmm-utils.R"))
		assignInNamespace("getNumberOfSNPs", getNumberOfSNPs.chip, "crlmm")

		XY = readGenCallOutput(file="chip.txt", cdfName=cdfName, colnames=list("SampleID"="Sample ID", "SNPID"="SNP Name", "XRaw"="X Raw", "YRaw"="Y Raw"), type=list("SampleID"="character", "SNPID"="character", "XRaw"="integer", "YRaw"="integer"))
	}
 
	crlmmOut <- try(crlmm::crlmmIllumina(XY=XY, verbose=FALSE, SNRMin=SNRMin))
	
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
	colnames(genotypes)<-paste("chip.", colnames(genotypes), sep="")

	if(length(grep("Probability", ls(assayData(crlmmOut)))) > 0) {
		flags<-assayData(crlmmOut)[["callProbability"]] #callProbability(crlmmOut)
		flags2<-flags
		flags2[flags<=quantile(flags, 0.95)]<-"P"
		flags2[flags>quantile(flags, 0.95)]<-"M"
		colnames(flags2)<-paste("flags.", colnames(flags), sep="")
		dat2 <- data.frame(genotypes, flags2)
	} else {
		dat2 <- genotypes
	}

	chiptype<-c(rep("Illumina-SNP", ncol(genotypes)))
	sample<-colnames(genotypes)
	group<-c(rep("", ncol(genotypes)))
	write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)
}

# Writing data to disk
write.table(dat2, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
