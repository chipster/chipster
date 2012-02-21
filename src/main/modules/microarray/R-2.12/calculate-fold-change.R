# TOOL calculate-fold-change.R: "Calculate fold change" (Calculates a geometric or arithmetic average of gene expression for replicate chips and then
# calculates a difference or ratio between the averages. The output fold change can be represented either in log2 or linear scale. Note that
# the tool is only applicable if you have exactly two groups of samples, and the reference group is determined by the order-parameter, not by a smaller group number.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT fold-change.tsv: fold-change.tsv 
# PARAMETER column: column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to average.)
# PARAMETER geometric TYPE [yes, no] DEFAULT yes (Should the geometric or arithmetic mean be used in the calculation of average expression for the sample groups?)
# PARAMETER scale TYPE [log2, linear] DEFAULT log2 (What scale to use for expressing the results. Log2 yields a symmetric distribution around zero with no change being equal to 0, up-regulation taking positive values and down-regulation negative values. Conversely, in linear scale 
# up-regulation is represented by values greater than 1 and down-regulation values being between 0 and 1.)
# PARAMETER order TYPE [normal, reversed] DEFAULT normal (Defines which samples to treat as the reference group. In the normal mode the group belonging
# of the first sample is assumed to be the reference, whereas in the reversed mode the situation is reversed. Note that this is different to how
# FC values are calculated in the various statistical testing tool, where the reference group is defined by the lower group number or alphabetical order.)

# Calculate fold changes between groups of samples
# JTT 30.7.2007
# MG, 3.5.2011, added parameters for choosing aritmetic or geometric mean
# and for choosing linear or log scale
# OH, 30.11.2011, support for names as group categories
# MG, 30.11.2011, added parameter do define reference group

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Test needs a parameter "groups" that specifies the grouping of the samples
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
groups<-phenodata[,pmatch(column,colnames(phenodata))]

# Sanity checks
if(length(unique(groups))==1) {
	stop("CHIPSTER-NOTE: You do not have any replicates to average!")
}
if(length(unique(groups))>2) {
	stop("CHIPSTER-NOTE: You have more than two groups! I don't know how to calculate fold change.")
}

# If arithmetic mean, then transform values to linear scale
if (geometric == "no") {
	dat2 <- as.data.frame (2^dat2)
}

# Calculating averages
columnnames<-c()
dat3<-matrix(nrow=nrow(dat2), ncol=length(unique(groups)), NA)
#for(i in 1:length(unique(groups))) {
index=1
for( i in unique(groups) ) {
	dat3[,index]<-rowSums(data.frame(dat2[,which(groups==i)]))/ncol(data.frame(dat2[,which(groups==i)]))
	columnnames<-c(columnnames, paste("group", i, sep=""))
	index=index+1
}
columnnames<-paste("chip.", columnnames, sep="")
colnames(dat3)<-columnnames
rownames(dat3)<-rownames(dat2)

# Put the reference group values in the first column
if (order=="reversed") dat3 <- dat3[,order(c(2,1))]

# Calculating the fold change
# Treatment divided by the control
if (geometric == "yes") {
	FC <- dat3[,2]-dat3[,1]
} else {
	FC <- dat3[,2] / dat3 [,1]
}

# If arithmetic mean and log2 scale, then transform values back to log2 scale
if (geometric == "no" && scale == "log2") {
	FC <- log2(FC)
}

# If geometric mean and linear scale, then transform values to linear scale
if (geometric == "yes" && scale == "linear") {
	FC <- 2^FC
}

# Saving the results
write.table(data.frame(dat, FC=FC), file="fold-change.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# EOF
