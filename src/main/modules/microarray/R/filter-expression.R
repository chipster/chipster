# TOOL filter-expression.R: "Filter by expression" (Filtering the genes by their expression. Note that chips are by default NOT scaled to have a mean of zero, which is most appropriate when applied to 1-color data, where this tool can be used for filtering out low quality or saturated probes. For 2-color data, the tool is most useful for filtering out invariant genes.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT expression-filter.tsv: expression-filter.tsv 
# PARAMETER over.expressed.cutoff: over.expressed.cutoff TYPE DECIMAL FROM -100000 TO 100000 DEFAULT 100 (Cut-off for over-expressed genes. For 1-color data this value should be a positive number, useful for excluding probes whoose expression level has reached the saturation limit. For 2-color data this threshold can be used to determine the lowest level of up-regulation.)
# PARAMETER under.expressed.cutoff: under.expressed.cutoff TYPE DECIMAL FROM -100000 TO 100000 DEFAULT 5 (Cut-off for under-expressed genes. For 1-color data this value should be a positive number, useful for excluding probes whoose expression level is under or near background level. For 2-color data this threshold can be used to determine the lowest level of down-regulation.)
# PARAMETER number.of.chips: number.of.chips TYPE INTEGER FROM 1 TO 10000 DEFAULT 1 (Number of chips)
# PARAMETER scale.to.same.mean: scale.to.same.mean TYPE [yes: yes, no: no] DEFAULT no (Scale the data to the same mean before filtering. Scaling to the mean is only recommended for 2-color arrays, where this filtering tool can be used to filter out invariant genes.)
# PARAMETER included.genes: included.genes TYPE [outside-the-range: outside-the-range, inside-the-range: inside-the-range] DEFAULT inside-the-range (Filtering method)


# Parameter settings (default) for testing purposes
#over.expressed.cutoff <-1
#under.expressed.cutoff <--1
#number.of.chips <-9
#scale.to.same.mean <-"yes"
#included.genes <-"outside-the-range"


# Filtering by expression
# JTT 9.6.2006
#
# modified by MG, 24.6.2010

# Loading the libraries
library(genefilter)

# Renaming variables
up<-over.expressed.cutoff
down<-under.expressed.cutoff
p<-number.of.chips
meth<-included.genes

# Reading data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Scaling the data to the same mean
if(scale.to.same.mean=="yes") {
   scaled.dat<-genescale(dat2)
} else {
   scaled.dat<-dat2
}

# Filter
# If a data value is smaller than up AND larger than down, recode it with 0, otherwise with 1
if(meth=="outside-the-range") {
	scaled.dat2<-ifelse(scaled.dat<=up & scaled.dat>=down, 0, 1)
}
# If a data value is smaller than up AND larger than down, recode it with 0, otherwise with 1
if(meth=="inside-the-range") {
	scaled.dat2<-ifelse(scaled.dat<=up & scaled.dat>=down, 1, 0)
}

# Calculate a sum of the number of chips the filtering criterion is fulfilled
s<-apply(as.data.frame(scaled.dat2), MARGIN=1, FUN="sum")

# Select only the rows for which the rowsum exceeds the number of chips (p)
#if(meth=="outside-the-range") {
#   dat2<-dat[which(s>=p),]
#}
#if(meth=="inside-the-range") {
#   dat2<-dat[which(s<p),] 
#}
dat2<-dat[which(s>=p),]

# Writing out the data
write.table(dat2, "expression-filter.tsv", sep="\t", row.names=T, col.names=T, quote=F)
