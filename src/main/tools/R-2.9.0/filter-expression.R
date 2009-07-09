# ANALYSIS Preprocessing/"Filter by expression" (Filtering the genes by their expression. 
# Note that chips are by default first scaled to have a mean of zero, regardless of the previous normalizations.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT expression-filter.tsv
# PARAMETER over.expressed.cutoff DECIMAL FROM -100000 TO 100000 DEFAULT 1 (Cut-off for over-expressed genes)
# PARAMETER under.expressed.cutoff DECIMAL FROM -100000 TO 100000 DEFAULT -1 (Cut-off for under-expressed genes)
# PARAMETER number.of.chips INTEGER FROM 1 TO 10000 DEFAULT 2 (Number of chips)
# PARAMETER scale.to.same.mean [yes, no] DEFAULT yes (Scale the data to the same mean before filtering)
# PARAMETER included.genes [outside-the-range, inside-the-range] DEFAULT outside-the-range (Filtering method)


#over.expressed.cutoff <-1
#under.expressed.cutoff <--1
#number.of.chips <-9
#scale.to.same.mean <-"yes"
#included.genes <-"outside-the-range"


# Filtering by expression
# JTT 9.6.2006

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
# If a data value is smaller than up or larger than down, recode it with 0, otherwise with 1
scaled.dat2<-ifelse(scaled.dat<=up & scaled.dat>=down, 0, 1)

# Calculate a sum of the number of chips the filtering criterion is fulfilled
s<-apply(as.data.frame(scaled.dat2), MARGIN=1, FUN="sum")

# Select only the rows for which the rowsum exceeds the number of chips (p)
if(meth=="outside-the-range") {
   dat2<-dat[which(s>=p),]
}
if(meth=="inside-the-range") {
   dat2<-dat[which(s<p),] 
}

# Writing out the data
write.table(dat2, "expression-filter.tsv", sep="\t", row.names=T, col.names=T, quote=F)
