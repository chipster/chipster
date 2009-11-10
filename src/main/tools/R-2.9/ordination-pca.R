# ANALYSIS Statistics/PCA (Principal component analysis. Creates principal components using genes or samples. 
# The number of principal component to save is controlled through the explained variablity. All principal components
# are saved until the explained variability is exceeded, but at least 3 components are always saved.
# Components are saved as new variables that can be used in further analyses.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT pca.tsv
# PARAMETER do.pca.on [genes, chips] DEFAULT genes (Data to use for the analysis)
# PARAMETER explained.variation PERCENT DEFAULT 80 (Amount of experimental variation to explain)
# PARAMETER scaling [yes, no] DEFAULT no (Scale the data to have a unit variance)
# PARAMETER centering [yes, no] DEFAULT yes (Scale the data to have the same mean)


# Principal component analysis
# Calculates principal components, and saves them for further analyses
# JTT 27.6.2006

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

# Saving the PCs with data
write.table(data.frame(round(pcs, digits=2)), "pca.tsv", sep="\t", row.names=T, col.names=T, quote=F)
