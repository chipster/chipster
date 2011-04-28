# TOOL generate-phenodata.R: "Generate phenodata" (If run on a prenormalized file, generates a blank phenodata for it.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# OUTPUT normalized.tsv: normalized.tsv 
# PARAMETER chiptype: chiptype TYPE STRING DEFAULT empty ()


# Combines two different tables using gene names
# JTT 6.7.2006
#
# modified by MG 19.1.2010

# Parameter settings (default) for testing purposes
#chiptype<-c("hgu133a")

# Figure out how the data is organized and load it
file<-c("normalized.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1, nrows=1)
ind.calls <- grep("flag", names(dat))
ind.flags <- grep("chip", names(dat))
ind.dat <- grep("chip", names(dat))
colclasses <- c(rep("numeric", length(ind.dat)), rep("character", 
				length(ind.calls)))
number.annotations <- length(names(dat))-length(ind.calls)-length(ind.flags)
colclasses <- append(rep("factor",1+number.annotations),colclasses)
dat <- read.table(file, header=T, sep="\t", row.names=1, colClasses=colclasses)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-data.frame(dat[,grep("chip", names(dat))])

# Check whether column names contain "chip." and, if so, remove from name
sample<-colnames(dat2)
if (length (grep("chip.", sample)) >= 1) {
	sample<-sub("chip.","",sample)
}

# Generates the variables
# sample<-colnames(dat2)
group<-c(rep("", ncol(dat2)))
if(chiptype=="empty") {
   chiptype<-c("empty")
} else {
   chiptype<-chiptype
}

# Writes out the data and the phenodata table
write.table(dat, file="normalized.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)
