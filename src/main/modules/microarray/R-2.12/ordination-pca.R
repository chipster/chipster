# TOOL ordination-pca.R: PCA (Principal component analysis. Creates principal components using genes or samples. The number of principal component to save is controlled through the explained variablity. All principal components are saved until the explained variability is exceeded, but at least 3 components are always saved. Components are saved as new variables that can be used in further analyses.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT pca.tsv: pca.tsv 
# PARAMETER do.pca.on: do.pca.on TYPE [genes: genes, chips: chips] DEFAULT genes (Data to use for the analysis)
# PARAMETER explained.variation: explained.variation TYPE PERCENT DEFAULT 80 (Amount of experimental variation to explain)
# PARAMETER scaling: scaling TYPE [yes: yes, no: no] DEFAULT no (Scale the data to have a unit variance)
# PARAMETER centering: centering TYPE [yes: yes, no: no] DEFAULT yes (Scale the data to have the same mean)


# Principal component analysis
# Calculates principal components, and saves them for further analyses
# JTT 27.6.2006
#
# modified by MG, 29.3.2010 to enable coloring by cluster feature 

# Parameter settings (default) for testing purposes
do.pca.on<-c("chips")
explained.variation<-c(80)
scaling<-c("no")
centering<-c("yes")

# Renaming variables
pcaon<-do.pca.on
expvar<-explained.variation

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# PCA calculations
if(scaling=="yes" & centering=="yes") {
	if(pcaon=="genes") {
		pc<-prcomp(dat2, scale=T, center=T)
	}
	if(pcaon=="chips") {
		# Transpose the matrix for sample-wise PCA 
		pc<-prcomp(t(dat2), scale=T, center=T)
	}
}
if(scaling=="yes" & centering=="no") {
	if(pcaon=="genes") {
		pc<-prcomp(dat2, scale=T, center=F)
	}
	if(pcaon=="chips") {
		# Transpose the matrix for sample-wise PCA 
		pc<-prcomp(t(dat2), scale=T, center=F)
	}
}
if(scaling=="no" & centering=="yes") {
	if(pcaon=="genes") {
		pc<-prcomp(dat2, scale=F, center=T)
	}
	if(pcaon=="chips") {
		# Transpose the matrix for sample-wise PCA 
		pc<-prcomp(t(dat2), scale=F, center=T)
	}
}
if(scaling=="no" & centering=="no") {
	if(pcaon=="genes") {
		pc<-prcomp(dat2, scale=F, center=F)
	}
	if(pcaon=="chips") {
		# Transpose the matrix for sample-wise PCA 
		pc<-prcomp(t(dat2), scale=F, center=F)
	}
}

# How many PCs to save?
no<-as.vector(head(which(summary(pc)$importance[3,]>=(expvar/100)), n=1)[1])
if(no<3) {
	no<-c(3)
}
if(pcaon=="genes") {
	pcs<-pc$x[,1:no]
}
if(pcaon=="chips") {
	pcs<-pc$x[,1:no]
}

# Converting PCs from matrix format into data frame
pcs<-as.data.frame(pcs)

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
