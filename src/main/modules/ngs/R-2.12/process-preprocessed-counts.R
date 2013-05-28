# TOOL process-preprocessed-counts.R: "Preprocess count table" (Converts count table to Chipster format and creates a phenodata file for it. You need to import your count table using the Import tool. Make sure to mark every count column as Sample and the column containing feature names as Identifier.)
# INPUT ngs{...}.tsv: ngs{...}.tsv TYPE CDNA 
# OUTPUT counts.tsv: counts.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv
# PARAMETER experimentype: "Experiment type" TYPE STRING DEFAULT empty (You can define experiment type, for example rna-seq.)
# PARAMETER OPTIONAL keep.annotations: "Keep annotations" TYPE [yes: yes, no: no] DEFAULT no (Keep or discard annotation columns after preprocessing.)


# MK 08.05.2013 Created based on the process prenormalized affy script


# Loads the libraries
library(limma)

# Reading data
columns<-list(R="sample")
annotation<-c("identifier")
columns.other<-c("annotation", "Annotation", "ANNOTATION")

files<-dir()
files<-files[files!="phenodata.tsv"]
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

# Mock normalization
dat2<-normalizeBetweenArrays(dat$E, method="none")
sample<-paste(colnames(dat2), ".tsv", sep="")
colnames(dat2)<-paste("chip.", colnames(dat2), sep="")

# Writes out a phenodata table
# sample<-paste(sample, ".tsv", sep="")
org_name <- files
group <- c(rep("", length(sample)))
chiptype <-"not applicable"
experiment <- c(rep(experimentype, length(sample)))
libsize <-c(rep("", length(sample)))

#sample, original name, chiptype, experiment, group, library size, description

write.table(data.frame(sample=sample, chiptype=chiptype, experiment=experiment, group=group, library_size=libsize, description=sample), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Writes the results into a file
if (keep.annotations=="yes") {
	output_table <- data.frame(symbol=annotation.info, dat2);	
} else {
	output_table <- data.frame(dat2);		
}
	
rownames (output_table) <- unlist(dat$genes)
write.table(output_table, file="counts.tsv", col.names=T, quote=F, sep="\t", row.names=T)
