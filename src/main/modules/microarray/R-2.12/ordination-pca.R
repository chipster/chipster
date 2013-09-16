# TOOL ordination-pca.R: PCA (Principal component analysis. Creates principal components using genes or samples. The number of principal component to save is controlled through the explained variablity. All principal components are saved until the explained variability is exceeded, but at least 3 components are always saved. Components are saved as new variables that can be used in further analyses.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT pca.tsv: pca.tsv 
# OUTPUT variance.tsv: variance.tsv 
# PARAMETER do.pca.on: do.pca.on TYPE [chips: chips] DEFAULT chips (Data to use for the analysis)
# PARAMETER explained.variation: explained.variation TYPE PERCENT DEFAULT 80 (Amount of experimental variation to explain)
# PARAMETER scaling: scaling TYPE [yes: yes, no: no] DEFAULT no (Scale the data to have a unit variance)
# PARAMETER centering: centering TYPE [yes: yes, no: no] DEFAULT yes (Scale the data to have the same mean)

# JTT, 27.6.2006
# MG, 29.3.2010 to enable coloring by cluster feature 
# IS, 12.10.2012 to cope with tables with gene descriptions (that typically contain 's)
# MK, 16.09.2013 modified to produce variance explanation table

# Parameter settings (default) for testing purposes
# do.pca.on<-c("chips")
# explained.variation<-c(80)
# scaling<-c("no")
# centering<-c("yes")

# Renaming variables
pcaon<-do.pca.on
expvar<-explained.variation

# Loads the normalized data
file <- 'normalized.tsv'
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, check.names=FALSE)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

if(pcaon=="chips") {
	dat2 <- t(dat2)
}

# PCA calculations
if(scaling=="yes" & centering=="yes") {
	pc<-prcomp(dat2, scale=T, center=T)
}
if(scaling=="yes" & centering=="no") {
	pc<-prcomp(dat2, scale=T, center=F)
}
if(scaling=="no" & centering=="yes") {
	pc<-prcomp(dat2, scale=F, center=T)
}
if(scaling=="no" & centering=="no") {
	pc<-prcomp(dat2, scale=F, center=F)
}

# How many PCs to save?
no<-as.vector(head(which(summary(pc)$importance[3,]>=(expvar/100)), n=1)[1])
if(no<3) {
	no<-c(3)
}
# Converting PCs from matrix format into data frame
pcs<-as.data.frame(pc$x[,1:no])

# Giving the PC headers new names
names(pcs)<-paste("chip.", names(pcs), sep="")
pcs<-round(pcs,digits=2)

# Add expression and other column values
if (pcaon=="genes") {
	dat4 <- cbind(dat,pcs)
}
if (pcaon=="chips") {
	dat4 <- pcs
}

# Saving the PCs with data
write.table(data.frame(dat4), "pca.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# Variance explained
write.table(round(data.frame(summary(pc)$importance[,1:no]),3), "variance.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# EOF
