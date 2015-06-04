# TOOL norm-prenormalized-affy.R: "Process prenormalized affy" (If you import prenormalized Affymetrix data that is not in Chipster format, you need to import it through the Import tool, and then use this tool for preprocessing it. During data import, make sure to mark every column containing normalized expression values as Sample and the column containing the Affymetrix probe ID:s as Identifier. If you want to be able to use annotation, you need to SPECIFY THE CHIPTYPE, e.g. hgu133a2.db.)
# INPUT microarray{...}.tsv: microarray{...}.tsv TYPE CDNA 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER chiptype: "Chiptype" TYPE STRING DEFAULT empty ()
# PARAMETER keep.annotations: "Keep annotations" TYPE [yes: yes, no: no] DEFAULT no (Keep or discard annotation column after preprocessing. Please note that gene symbol information associated with the given chiptype will be replaced by data-specific annotations by setting this parameter to "yes")
# PARAMETER keep.flags: "Keep flags" TYPE [yes: yes, no: no] DEFAULT no (Keep or discard flag-columns after preprocessing. Please note that flag columns must have been named as "flag")

# Process prenormalized affy
# MG 6.4.2010
# MK 18.2.2013
# EK 23.4.2013 added .tsv ending to expression column names so that sample renaming is possible in interactive visualizations

# Loads the libraries
library(limma)
library(affy)
library(gcrma)

# Reading data
columns<-list(R="sample")
annotation<-c("identifier")
columns.other<-c("flag", "Flag", "FLAG", "annotation", "Annotation", "ANNOTATION")

files<-dir(pattern = "microarray")
#seems to check columns from the first file only. If the column is not there, attribute is ignored. If the column is there, but not
#in the others, function crashes
dat<-read.maimages(files=files, columns=columns, annotation=annotation, other.columns=columns.other) 

#dat <- read.table(file="normalized.tsvt", header=TRUE, sep="\t")

# Extract annotations
annotation.info <- NULL;
if (keep.annotations=="yes") {
	if(length(grep("^annotation$", unlist(attributes(dat$other)), ignore.case=T))) {
		annotation.col <- unlist(attributes(dat$other))[grep("^annotation$", unlist(attributes(dat$other)), ignore.case=T)]
		annotation.info	<- dat$other[[annotation.col]];
		annotation.info <- as.character(annotation.info[,1]);
	}
}

flag.info <- NULL;
if (keep.flags=="yes") {
	if(length(grep("^flag$", unlist(attributes(dat$other)), ignore.case=T))) {
		flag.col 			<- unlist(attributes(dat$other))[grep("^flag$", unlist(attributes(dat$other)), ignore.case=T)]
		flag.info 			<- dat$other$flag;
		colnames(flag.info)	<-paste("flag.", colnames(flag.info), sep="")
	}
}

# Mock normalization
dat2<-normalizeBetweenArrays(dat$E, method="none")
sample<-colnames(dat2)
colnames(dat2)<-paste("chip.", colnames(dat2), sep="")

# Writes out a phenodata table
# sample<-paste(sample, ".tsv", sep="")
group<-c(rep("", length(sample)))
training<-c(rep("", length(sample)))
time<-c(rep("", length(sample)))
random<-c(rep("", length(sample)))
if(chiptype=="empty") {
	chiptype<-"Affy"
}
write.table(data.frame(sample=sample, chiptype=chiptype, group=group, description=sample), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Writing out data
a<-try(library(chiptype, character.only=T))
if(class(a)=="try-error") { a<-try(library(paste(chiptype, ".db", sep=""), character.only=T)) }

if(chiptype!="empty" & class(a)!="try-error") {
	# Including gene names to data
	lib2		<- sub('.db','',chiptype)
	symbol		<- gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep="")))))[unlist(dat$genes),])
	genename	<- gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep="")))))[unlist(dat$genes),])
	genename 	<- gsub("#", "", genename)
	
	# Writes the results into a file
	if (keep.flags=="yes" & keep.annotations=="yes") {
		output_table <- data.frame(symbol=annotation.info, description=genename, dat2, flag.info);
	} else if (keep.flags=="yes" & keep.annotations=="no") {
		output_table <- data.frame(symbol, description=genename, dat2, flag.info);
	} else if (keep.flags=="no" & keep.annotations=="yes") {
		output_table <- data.frame(symbol=annotation.info, description=genename, dat2);
	} else {
		output_table <- data.frame(symbol, description=genename, dat2);		
	}
	
	rownames (output_table) <- make.names(unlist(dat$genes), unique=T)
	rownames(output_table) <- gsub("^X", "", rownames(output_table))
	write.table(output_table, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
}

if(chiptype=="empty" | class(a)=="try-error") {
	
	# Writes the results into a file
	if (keep.flags=="yes" & keep.annotations=="yes") {
		output_table <- data.frame(symbol=annotation.info, dat2, flag.info);
	} else if (keep.flags=="yes" & keep.annotations=="no") {
		output_table <- data.frame(dat2, flag.info);
	} else if (keep.flags=="no" & keep.annotations=="yes") {
		output_table <- data.frame(symbol=annotation.info, dat2);	
	} else {
		output_table <- data.frame(dat2);		
	}
	
	rownames (output_table) <- make.names(unlist(dat$genes), unique=T)
	rownames(output_table) <- gsub("^X", "", rownames(output_table))
	write.table(output_table, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
}

