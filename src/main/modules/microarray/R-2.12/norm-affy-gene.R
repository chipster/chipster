# TOOL norm-affy-gene.R: "Affymetrix gene arrays" (Affymetrix RMA preprocessing for CEL-files. Please note that the preprocessing might take a long time. YOU HAVE TO SPECIFY THE CHIPTYPE.)
# INPUT microarray{...}.cel: microarray{...}.cel TYPE AFFY 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER chiptype: chiptype TYPE [empty: empty, human-1.0-ST: human-1.0-ST, human-1.1-ST: human-1.1-ST, human-2.0-ST: human-2.0-ST, human-2.1-ST: human-2.1-ST, mouse-1.0-ST: mouse-1.0-ST, mouse-1.1-ST: mouse-1.1-ST, mouse-2.0-ST: mouse-2.0-ST, mouse-2.1-ST: mouse-2.1-ST, rat-1.0-ST: rat-1.0-ST, rat-1.1-ST: rat-1.1-ST, rat-2.0-ST: rat-2.0-ST, rat-2.1-ST: rat-2.1-ST, zebra_fish-1.0-ST: zebra_fish-1.0-ST, zebra_fish-1.1-ST: zebra_fish-1.1-ST, arabidopsis-1.0-ST-entrez: arabidopsis-1.0-ST-entrez, arabidopsis-1.1-ST-entrez: arabidopsis-1.1-ST-entrez, arabidopsis-1.0-ST-tair: arabidopsis-1.0-ST-tair, arabidopsis-1.1-ST-tair: arabidopsis-1.1-ST-tair] DEFAULT empty (Chiptype)
# PARAMETER biomart: "biomaRt annotation" TYPE [yes: annotate, no: skip] DEFAULT no (In the case where no annotation has been attached to CDF-files, attach symbol and description information to probesets using bioMart)

# Affymetrix normalization
# JTT, 3.2.2009
# MG. 18.10.2011, added version 1.1 ST arrays
# MK 15.05.2014 added new arrays

# Initializes analyses
library(affy)
library(biomaRt)

# Reads in data
dat<-ReadAffy()
if(chiptype=="empty") {
	stop("CHIPSTER-NOTE: You need to specify the chiptype. Please run the tool again.")
}

filt <- "entrezgene"
usemart <- "ensembl"

if(chiptype=="human-1.0-ST") {
	dat@cdfName<-"hugene10stv1hsentrezgcdf"
	dat@annotation<-"hugene10stv1hsentrezgcdf"
	chiptype<-"hugene10stv1hsentrezg.db"
	dataset <- "hsapiens_gene_ensembl"
}
if(chiptype=="human-1.1-ST") {
	dat@cdfName<-"hugene11stv1hsentrezgcdf"
	dat@annotation<-"hugene11stv1hsentrezgcdf"
	chiptype<-"hugene11stv1hsentrezg.db"
	dataset <- "hsapiens_gene_ensembl"
}
if(chiptype=="human-2.0-ST") {
	dat@cdfName<-"hugene20sthsentrezgcdf"
	dat@annotation<-"hugene20sthsentrezgcdf"
	chiptype<-"hugene20sthsentrezg.db"
	dataset <- "hsapiens_gene_ensembl"
}
if(chiptype=="human-2.1-ST") {
	dat@cdfName<-"hugene21sthsentrezgcdf"
	dat@annotation<-"hugene21sthsentrezgcdf"
	chiptype<-"hugene21sthsentrezg.db"
	dataset <- "hsapiens_gene_ensembl"
}
if(chiptype=="mouse-1.0-ST") {
	dat@cdfName<-"mogene10stv1mmentrezgcdf"
	dat@annotation<-"mogene10stv1mmentrezgcdf"
	chiptype<-"mogene10stv1mmentrezg.db"
	dataset <- "mmusculus_gene_ensembl"
}
if(chiptype=="mouse-1.1-ST") {
	dat@cdfName<-"mogene11stv1mmentrezgcdf"
	dat@annotation<-"mogene11stv1mmentrezgcdf"
	chiptype<-"mogene11stv1mmentrezg.db"
	dataset <- "mmusculus_gene_ensembl"
}
if(chiptype=="mouse-2.0-ST") {
	dat@cdfName<-"mogene20stmmentrezgcdf"
	dat@annotation<-"mogene20stmmentrezgcdf"
	chiptype<-"mogene20stmmentrezg.db"
	dataset <- "mmusculus_gene_ensembl"
}
if(chiptype=="mouse-2.1-ST") {
	dat@cdfName<-"mogene21stmmentrezgcdf"
	dat@annotation<-"mogene21stmmentrezgcdf"
	chiptype<-"mogene21stmmentrezg.db"
	dataset <- "mmusculus_gene_ensembl"
}
if(chiptype=="rat-1.0-ST") {
	dat@cdfName<-"ragene10stv1rnentrezgcdf"
	dat@annotation<-"ragene10stv1rnentrezgcdf"
	chiptype<-"ragene10stv1rnentrezg.db"
	dataset <- "rnorvegicus_gene_ensembl"
}
if(chiptype=="rat-1.1-ST") {
	dat@cdfName<-"ragene11stv1rnentrezgcdf"
	dat@annotation<-"ragene11stv1rnentrezgcdf"
	chiptype<-"ragene11stv1rnentrezg.db"
	dataset <- "rnorvegicus_gene_ensembl"
}
if(chiptype=="rat-2.0-ST") {
	dat@cdfName<-"ragene20strnentrezgcdf"
	dat@annotation<-"ragene20strnentrezgcdf"
	chiptype<-"ragene20strnentrezg.db"
	dataset <- "rnorvegicus_gene_ensembl"
}
if(chiptype=="rat-2.1-ST") {
	dat@cdfName<-"ragene21strnentrezgcdf"
	dat@annotation<-"ragene21strnentrezgcdf"
	chiptype<-"ragene21strnentrezg.db"
	dataset <- "rnorvegicus_gene_ensembl"
}
if(chiptype=="zebra_fish-1.0-ST") {
	dat@cdfName<-"zebgene10stdrentrezgcdf"
	dat@annotation<-"zebgene10stdrentrezgcdf"
	chiptype<-"zebgene10stdrentrezg.db"
	dataset <- "drerio_gene_ensembl"
}
if(chiptype=="zebra_fish-1.1-ST") {
	dat@cdfName<-"zebgene11stdrentrezgcdf"
	dat@annotation<-"zebgene11stdrentrezgcdf"
	chiptype<-"zebgene11stdrentrezg.db"
	dataset <- "drerio_gene_ensembl"
}
if(chiptype=="arabidopsis-1.0-ST-entrez") {
	dat@cdfName<-"aragene10statentrezgcdf"
	dat@annotation<-"aragene10statentrezgcdf"
	chiptype<-"aragene10statentrezg.db"
	usemart <- as.character(listMarts()[grep("plants_mart", listMarts()[,1]),1])
	dataset <- "athaliana_eg_gene"
}
if(chiptype=="arabidopsis-1.1-ST-entrez") {
	dat@cdfName<-"aragene11statentrezgcdf"
	dat@annotation<-"aragene11statentrezgcdf"
	chiptype<-"aragene11statentrezg.db"
	usemart <- as.character(listMarts()[grep("plants_mart", listMarts()[,1]),1])
	dataset <- "athaliana_eg_gene"
}
if(chiptype=="arabidopsis-1.0-ST-tair") {
	dat@cdfName<-"aragene10stattairgcdf"
	dat@annotation<-"aragene10stattairgcdf"
	chiptype<-"aragene10stattairg.db"
	usemart <- as.character(listMarts()[grep("plants_mart", listMarts()[,1]),1])
	dataset <- "athaliana_eg_gene"
	filt <- "ensembl_gene_id"
}
if(chiptype=="arabidopsis-1.1-ST-tair") {
	dat@cdfName<-"aragene11stattairgcdf"
	dat@annotation<-"aragene11stattairgcdf"
	chiptype<-"aragene11stattairg.db"
	usemart <- as.character(listMarts()[grep("plants_mart", listMarts()[,1]),1])
	dataset <- "athaliana_eg_gene"
	filt <- "ensembl_gene_id"
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
		
		for (i in 1:ncol(annotated_dataf)) {
			annotated_dataf[,i] <- gsub("\'+", "", annotated_dataf[,i])
			annotated_dataf[,i] <- gsub("\"+", "", annotated_dataf[,i])
			annotated_dataf[,i] <- gsub("\\#+", "", annotated_dataf[,i])
		}

		write.table(data.frame(symbol=annotated_dataf[gene_id, "external_gene_id"], description=annotated_dataf[gene_id, "description"], dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	} else {
		write.table(data.frame(dat2), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	}
}


