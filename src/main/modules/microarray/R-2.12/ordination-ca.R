# TOOL ordination-ca.R: DCA (Does a detrended correspondence analysis for the data. Can be used for, e.g., quality control.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT ca.png: ca.png 
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Non-metric multidimensional scaling
# JTT 30.1.2009

# Parameter settings (default) for testing purposes
#image.width<-c(600)
#image.height<-c(600)

# Renaming variables
w<-image.width
h<-image.height

# Loads the libraries
library(vegan)
library(MASS)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Decorana can't handle negative values
if(any(dat2<0)) {
   dat2<-na.omit(dat2)
   dat2<-dat2+abs(min(dat2))
}

# Detrended correspondence analysis
fit<-decorana(dat2)

# Plotting the image
bitmap(file="ca.png", width=w/72, height=h/72)
plot(fit, main="Detrended correspondence analysis", type="n")
points(fit, display=c("species"), pch=15, col=phenodata$group)
points(fit, display=c("sites"), pch=20, col=c("black"), cex=0.6)
text(fit, display=c("species"), labels=phenodata$description, cex=0.75, adj=1)
dev.off()
