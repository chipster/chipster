# ANALYSIS Statistics/NMDS (Non-metric multidimensional scaling. Creates a 2-dimensional representation of 
# the arrays using Euclidean distances. Can be used for quality control.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT nmds.pdf
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Non-metric multidimensional scaling
# JTT 14.1.2008

# MG, 15.9.2011
# Updated colors and legend

# Parameter settings (default) for testing purposes
#image.width<-c(600)
#image.height<-c(600)

# Renaming variables
w<-image.width
h<-image.height

# Loads the libraries
library(MASS)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Calculating the distance matrix
dat2.dist<-dist(t(dat2))

# Calculating the MDS
mds<-isoMDS(dat2.dist)

# Plotting the image
pdf(file="nmds.pdf", width=w/72, height=h/72)
plot(mds$points[,1], mds$points[,2], main="NMDS", pch=19, xlab="Dimension 1", ylab="Dimension 2", type="n")
text(mds$points[,1], mds$points[,2], phenodata$description, cex=0.75, col=phenodata$group+2)
if (length(levels(as.factor(phenodata$group))) > 0) {
	legend (x="topleft", legend=levels(as.factor(phenodata$group)), col=1+1:length(levels(as.factor(phenodata$group))), cex=0.5, pch=19)
}
dev.off()
