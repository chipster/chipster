# TOOL norm-affy-gene.R: "Affymetrix gene arrays" (Affymetrix RMA preprocessing for CEL-files. Please note that the preprocessing might take a long time. YOU HAVE TO SPECIFY THE CHIPTYPE.)
# INPUT microarray{...}.cel: microarray{...}.cel TYPE AFFY 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER chiptype: "Chiptype" TYPE [empty: empty, oligo: oligo, human-1.0-ST: human-1.0-ST, human-1.1-ST: human-1.1-ST, human-2.0-ST: human-2.0-ST, human-2.1-ST: human-2.1-ST, human-hta20: human-hta20, human-prime: human-prime, mouse-1.0-ST: mouse-1.0-ST, mouse-1.1-ST: mouse-1.1-ST, mouse-2.0-ST: mouse-2.0-ST, mouse-2.1-ST: mouse-2.1-ST, rat-1.0-ST: rat-1.0-ST, rat-1.1-ST: rat-1.1-ST, rat-2.0-ST: rat-2.0-ST, rat-2.1-ST: rat-2.1-ST, zebra_fish-1.0-ST: zebra_fish-1.0-ST, zebra_fish-1.1-ST: zebra_fish-1.1-ST, arabidopsis-1.0-ST-entrez: arabidopsis-1.0-ST-entrez, arabidopsis-1.1-ST-entrez: arabidopsis-1.1-ST-entrez, arabidopsis-1.0-ST-tair: arabidopsis-1.0-ST-tair, arabidopsis-1.1-ST-tair: arabidopsis-1.1-ST-tair, fly-1.0-ST: fly-1.0-ST, fly-1.1-ST: fly-1.1-ST, celegans-1.0-ST: celegans-1.0-ST, celegans-1.1-ST: celegans-1.1-ST, dog-1.0-ST: dog-1.0-ST, dog-1.1-ST: dog-1.1-ST, rice-1.0-ST-entrez: rice-1.0-ST-entrez, rice-1.1-ST-entrez: rice-1.1-ST-entrez] DEFAULT empty (Chiptype. The oligo class has platform design information for various Affymetrix chips, but at the moment corresponding annotation packages are available only for human, mouse and rat chips.)
# PARAMETER biomart: "biomaRt annotation" TYPE [yes: annotate, no: skip] DEFAULT no (In the case where no annotation has been attached to CDF-files, attach symbol and description information to probesets using bioMart)

# Affymetrix normalization
# JTT, 3.2.2009
# MG 18.10.2011, added version 1.1 ST arrays
# MK 15.05.2014, added new arrays
# MK 12.02.2014, added hta-20 chip
# MK 13.02.2014, annotation for oligo packages
# MK 28.03.2014, added dog, fly, c.elegans, human-prime chips

# Initializes analyses
library(oligo)
library(affy)
library(biomaRt)

# Reads in data
if(chiptype=="empty") {
	stop("CHIPSTER-NOTE: You need to specify the chiptype. Please run the tool again.")
}

filt <- "entrezgene"
usemart <- "ensembl"

if(chiptype=="human-prime") {
	custom_cdf = "primeviewhsentrezgcdf"
	chiptype<-"primeviewhsentrezg.db"
	dataset <- "hsapiens_gene_ensembl"
}
if(chiptype=="human-1.0-ST") {
	custom_cdf = "hugene10sthsentrezgcdf"
	#dat@cdfName<-"hugene10sthsentrezgcdf"
	#dat@annotation<-"hugene10sthsentrezgcdf"
	chiptype<-"hugene10sthsentrezg.db"
	dataset <- "hsapiens_gene_ensembl"
}
if(chiptype=="human-1.1-ST") {
	custom_cdf = "hugene11sthsentrezgcdf"
	#dat@cdfName<-"hugene11sthsentrezgcdf"
	#dat@annotation<-"hugene11sthsentrezgcdf"
	chiptype<-"hugene11sthsentrezg.db"
	dataset <- "hsapiens_gene_ensembl"
}
if(chiptype=="human-2.0-ST") {
	custom_cdf = "hugene20sthsentrezgcdf"
	#dat@cdfName<-"hugene20sthsentrezgcdf"
	#dat@annotation<-"hugene20sthsentrezgcdf"
	chiptype<-"hugene20sthsentrezg.db"
	dataset <- "hsapiens_gene_ensembl"
}
if(chiptype=="human-2.1-ST") {
	custom_cdf = "hugene21sthsentrezgcdf"
	#dat@cdfName<-"hugene21sthsentrezgcdf"
	#dat@annotation<-"hugene21sthsentrezgcdf"
	chiptype<-"hugene21sthsentrezg.db"
	dataset <- "hsapiens_gene_ensembl"
}
if(chiptype=="human-hta20") {
	custom_cdf = "hta20hsentrezgcdf"
	#dat@cdfName<-"hugene21sthsentrezgcdf"
	#dat@annotation<-"hugene21sthsentrezgcdf"
	chiptype<-"hta20hsentrezg.db"
	dataset <- "hsapiens_gene_ensembl"
}
if(chiptype=="mouse-1.0-ST") {
	custom_cdf = "mogene10stmmentrezgcdf"
	#dat@cdfName<-"mogene10stmmentrezgcdf"
	#dat@annotation<-"mogene10stmmentrezgcdf"
	chiptype<-"mogene10stmmentrezg.db"
	dataset <- "mmusculus_gene_ensembl"
}
if(chiptype=="mouse-1.1-ST") {
	custom_cdf = "mogene11stmmentrezgcdf"
	#dat@cdfName<-"mogene11stmmentrezgcdf"
	#dat@annotation<-"mogene11stmmentrezgcdf"
	chiptype<-"mogene11stmmentrezg.db"
	dataset <- "mmusculus_gene_ensembl"
}
if(chiptype=="mouse-2.0-ST") {
	custom_cdf = "mogene20stmmentrezgcdf"
	#dat@cdfName<-"mogene20stmmentrezgcdf"
	#dat@annotation<-"mogene20stmmentrezgcdf"
	chiptype<-"mogene20stmmentrezg.db"
	dataset <- "mmusculus_gene_ensembl"
}
if(chiptype=="mouse-2.1-ST") {
	custom_cdf = "mogene21stmmentrezgcdf"
	#dat@cdfName<-"mogene21stmmentrezgcdf"
	#dat@annotation<-"mogene21stmmentrezgcdf"
	chiptype<-"mogene21stmmentrezg.db"
	dataset <- "mmusculus_gene_ensembl"
}
if(chiptype=="rat-1.0-ST") {
	custom_cdf = "ragene10strnentrezgcdf"
	#dat@cdfName<-"ragene10strnentrezgcdf"
	#dat@annotation<-"ragene10strnentrezgcdf"
	chiptype<-"ragene10strnentrezg.db"
	dataset <- "rnorvegicus_gene_ensembl"
}
if(chiptype=="rat-1.1-ST") {
	custom_cdf = "ragene11strnentrezgcdf"
	#dat@cdfName<-"ragene11strnentrezgcdf"
	#dat@annotation<-"ragene11strnentrezgcdf"
	chiptype<-"ragene11strnentrezg.db"
	dataset <- "rnorvegicus_gene_ensembl"
}
if(chiptype=="rat-2.0-ST") {
	custom_cdf = "ragene20strnentrezgcdf"
	#dat@cdfName<-"ragene20strnentrezgcdf"
	#dat@annotation<-"ragene20strnentrezgcdf"
	chiptype<-"ragene20strnentrezg.db"
	dataset <- "rnorvegicus_gene_ensembl"
}
if(chiptype=="rat-2.1-ST") {
	custom_cdf = "ragene21strnentrezgcdf"
	#dat@cdfName<-"ragene21strnentrezgcdf"
	#dat@annotation<-"ragene21strnentrezgcdf"
	chiptype<-"ragene21strnentrezg.db"
	dataset <- "rnorvegicus_gene_ensembl"
}
if(chiptype=="zebra_fish-1.0-ST") {
	custom_cdf = "zebgene10stdrentrezgcdf"
	#dat@cdfName<-"zebgene10stdrentrezgcdf"
	#dat@annotation<-"zebgene10stdrentrezgcdf"
	chiptype<-"zebgene10stdrentrezg.db"
	dataset <- "drerio_gene_ensembl"
}
if(chiptype=="zebra_fish-1.1-ST") {
	custom_cdf = "zebgene11stdrentrezgcdf"
	#dat@cdfName<-"zebgene11stdrentrezgcdf"
	#dat@annotation<-"zebgene11stdrentrezgcdf"
	chiptype<-"zebgene11stdrentrezg.db"
	dataset <- "drerio_gene_ensembl"
}

if(chiptype=="fly-1.0-ST") {
	custom_cdf = "drogene10stdmentrezgcdf"
	chiptype<-"drogene10stdmentrezg.db"
	dataset <- "dmelanogaster_gene_ensembl"
}
if(chiptype=="fly-1.1-ST") {
	custom_cdf = "drogene11stdmentrezgcdf"
	chiptype<-"drogene11stdmentrezg.db"
	dataset <- "dmelanogaster_gene_ensembl"
}
if(chiptype=="celegans-1.0-ST") {
	custom_cdf = "elegene10stceentrezgcdf"
	chiptype<-"elegene10stceentrezg.db"
	dataset <- "celegans_gene_ensembl"
}
if(chiptype=="celegans-1.1-ST") {
	custom_cdf = "elegene11stceentrezgcdf"
	chiptype<-"elegene11stceentrezg.db"
	dataset <- "celegans_gene_ensembl"
}
if(chiptype=="dog-1.0-ST") {
	custom_cdf = "cangene10stcfentrezgcdf"
	chiptype<-"cangene10stcfentrezg.db"
	dataset <- "cfamiliaris_gene_ensembl"
}
if(chiptype=="dog-1.1-ST") {
	custom_cdf = "cangene11stcfentrezgcdf"
	chiptype<-"cangene11stcfentrezg.db"
	dataset <- "cfamiliaris_gene_ensembl"
}


if(chiptype=="arabidopsis-1.0-ST-entrez") {
	custom_cdf = "aragene10statentrezgcdf"
	#dat@cdfName<-"aragene10statentrezgcdf"
	#dat@annotation<-"aragene10statentrezgcdf"
	chiptype<-"aragene10statentrezg.db"
	usemart <- as.character(listMarts()[grep("plants_mart", listMarts()[,1]),1])
	#usemart <- as.character(listMarts()[grep("GRAMENE.*GENES", listMarts()[,2]),1])
	dataset <- "athaliana_eg_gene"
}
if(chiptype=="arabidopsis-1.1-ST-entrez") {
	custom_cdf = "aragene11statentrezgcdf"
	#dat@cdfName<-"aragene11statentrezgcdf"
	#dat@annotation<-"aragene11statentrezgcdf"
	chiptype<-"aragene11statentrezg.db"
	usemart <- as.character(listMarts()[grep("plants_mart", listMarts()[,1]),1])
	#usemart <- as.character(listMarts()[grep("GRAMENE.*GENES", listMarts()[,2]),1])
	dataset <- "athaliana_eg_gene"
}
if(chiptype=="arabidopsis-1.0-ST-tair") {
	custom_cdf = "aragene10stattairgcdf"
	#dat@cdfName<-"aragene10stattairgcdf"
	#dat@annotation<-"aragene10stattairgcdf"
	chiptype<-"aragene10stattairg.db"
	usemart <- as.character(listMarts()[grep("plants_mart", listMarts()[,1]),1])
	#usemart <- as.character(listMarts()[grep("GRAMENE.*GENES", listMarts()[,2]),1])
	dataset <- "athaliana_eg_gene"
	filt <- "ensembl_gene_id"
}
if(chiptype=="arabidopsis-1.1-ST-tair") {
	custom_cdf = "aragene11stattairgcdf"
	#dat@cdfName<-"aragene11stattairgcdf"
	#dat@annotation<-"aragene11stattairgcdf"
	chiptype<-"aragene11stattairg.db"
	usemart <- as.character(listMarts()[grep("plants_mart", listMarts()[,1]),1])
	#usemart <- as.character(listMarts()[grep("GRAMENE.*GENES", listMarts()[,2]),1])
	dataset <- "athaliana_eg_gene"
	filt <- "ensembl_gene_id"
}

if(chiptype=="rice-1.0-ST-entrez") {
	custom_cdf = "ricegene10stosentrezgcdf"
	usemart <- as.character(listMarts()[grep("plants_mart", listMarts()[,1]),1])
	dataset <- "osativa_eg_gene"
	filt <- "entrezgene"
}
if(chiptype=="rice-1.1-ST-entrez") {
	custom_cdf = "ricegene11stosentrezgcdf"
	usemart <- as.character(listMarts()[grep("plants_mart", listMarts()[,1]),1])
	dataset <- "osativa_eg_gene"
	filt <- "entrezgene"
}


# Normalizations
if(chiptype != "oligo") {
	dat2 <- justRMA(filenames=list.celfiles(), cdfname=custom_cdf)
} else {
	data.raw <- read.celfiles(filenames=list.celfiles())
	dat2 <- oligo::rma(data.raw)
}

dat2<-as.data.frame(round(exprs(dat2), digits=2))
names(dat2)<-paste("chip.", names(dat2), sep="")

if(chiptype == "oligo") {
	chiptype <- annotation(data.raw)
	chiptype <- gsub("pd\\.", "", chiptype)
	chiptype <- gsub("\\.st\\..*$", "", chiptype)
	chiptype <- gsub("\\.", "", chiptype)
	chiptype <- paste(chiptype, "sttranscriptcluster.db", sep="")
}

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
	# Including gene names to data
	lib2<-sub('.db','',chiptype)
	symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep="")))))[rownames(dat2),])
	genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep="")))))[rownames(dat2),])
	# Fxes an issue introduced in BioC2.4 where the "#" character is introduced in some gene names
	genename <- gsub("#", "", genename)
	symbol <- gsub("'", "", symbol)
	genename <- gsub("'", "", genename)

	# remove odd characters
	symbol <- gsub("\'|#|\"|\n|\t", "", symbol, perl=T) 
	genename <- gsub("\'|#|\"|\n|\t", "", genename, perl=T) 
	
	# Writes the results into a file
	write.table(data.frame(symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
} 

#
#In the case where no information is attached to the cdf, BioMart is used
if(chiptype=="empty" | class(a)=="try-error") {
	gene_id <- gsub("_at", "", rownames(dat2))
	if(length(gene_id) > 0 && biomart == "yes") {
		biom <- useMart(usemart, dataset=dataset)
		annotated_genes <- getBM(mart=biom, attributes=c(filt, "ensembl_gene_id", "external_gene_id","description"), filters=filt, values=gene_id,  uniqueRows = TRUE)
		annotated_genes <- annotated_genes[!duplicated(annotated_genes[,1]), ]
		
		annotated_dataf <- data.frame(annotated_genes)
		rownames(annotated_dataf) <- annotated_genes[, 1]

		# remove odd characters
		symbol <- gsub("\'|#|\"|\n|\t", "", annotated_dataf[gene_id, "external_gene_id"], perl=T) 
		genename <- gsub("\'|#|\"|\n|\t", "", annotated_dataf[gene_id, "description"], perl=T) 
		
		write.table(data.frame(symbol=symbol, description=genename, dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	} else {
		write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	}
}
