# TOOL process-preprocessed-counts.R: "Preprocess count table" (Converts count table to Chipster format and creates a phenodata file for it. You need to import your count table using the Import tool. Make sure to mark every count column as Sample and the column containing feature names as Identifier.)
# INPUT ngs{...}.tsv: ngs{...}.tsv TYPE CDNA 
# OUTPUT counts.tsv: counts.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv
# PARAMETER experimentype: "Experiment type" TYPE STRING DEFAULT empty (You can define experiment type, for example rna-seq.)
# PARAMETER OPTIONAL keep.annotations: "Keep annotations" TYPE [yes: yes, no: no] DEFAULT no (Keep or discard annotation columns after preprocessing.)

# MK 08.05.2013 Created based on the process prenormalized affy script
# MK 11.04.2013 Support for multiple annotation columns added

# Read files
#files<-dir()
#files<-files[files!="phenodata.tsv"]
files<-dir(pattern = "ngs")

# Loads the limma library
library(limma)

# Define accepted column headers
columns <- list(R="sample")
annotation <- c("identifier")
columns.other <- c("annotation", "Annotation", "ANNOTATION")

# Check how many annotation columns, fix name of annotation columns >2 and modify files accordingly
header <- t(read.table(files[1], header=F, nrows = 1))
anno.cols <- grep("annotation", t(header), ignore.case=T)
if(length(anno.cols) > 1) {
	header[anno.cols[-1]] <- paste(header[anno.cols[-1]], 1:length(anno.cols[-1]), sep="")
	replace.header <- paste(header, collapse = "\t")
	system(paste("find n*.tsv -exec sed -i -e '1s/.*/", replace.header,"/' {} \\; ", sep=""))
	columns.other <- c(columns.other, paste(columns.other, rep(1:length(anno.cols[-1]), each=3), sep=""))
}

# Read data
dat<-read.maimages(files=files, columns=columns, annotation=annotation, other.columns=columns.other) 

# Extract annotations
annotation.info <- NULL;
if (keep.annotations=="yes") {
	if(length(grep("^annotation", unlist(attributes(dat$other)), ignore.case=T)) > 0) {
		annotation.col <- unlist(attributes(dat$other))[grep("^annotation", unlist(attributes(dat$other)), ignore.case=T)]
		annotation.info	<- NULL
		for(i in 1:length(annotation.col)) {
			annotation.info	<- cbind(annotation.info, (dat$other[[annotation.col[i]]])[,1]);
		}
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
if (keep.annotations=="yes" && (!is.null(annotation.info))) {
	if(ncol(annotation.info) == 1) {
		output_table <- data.frame(symbol=annotation.info, dat2);	
	} else {
		colnames(annotation.info) <- c("Annotation", paste("Annotation", 2:(length(anno.cols[-1])+1), sep=""))
		output_table <- data.frame(annotation.info, dat2);
	}
} else {
	output_table <- data.frame(dat2);		
}
	
rownames (output_table) <- unlist(dat$genes)
write.table(output_table, file="counts.tsv", col.names=T, quote=F, sep="\t", row.names=T)
