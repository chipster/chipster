# ANALYSIS Visualisation/"Heatmap" (Draws a heatmap using Pearson correlation and average linkage.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT heatmap.png
# PARAMETER coloring.scheme [Green-Red, Blue-Yellow, Black-White] DEFAULT Blue-Yellow (Coloring scheme for the SOM map)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Heatmap
# JTT 3.10.2007

# Renaming variables
colpar<-coloring.scheme
w<-image.width
h<-image.height

# Loading the libraries
library(amap)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]
colnames(dat2)<-gsub(" ", "", phenodata$description)

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Does the clustering
clustg<-as.dendrogram(hcluster(x=dat2, method="pearson", link="average"))
clustc<-as.dendrogram(hcluster(x=t(dat2), method="pearson", link="average"))

# Generating the colors
if(colpar=="Green-Red") {
   heatcol<-colorRampPalette(c("Green", "Red"))(32)
}
if(colpar=="Blue-Yellow") {
   heatcol<-colorRampPalette(c("Blue", "Yellow"))(32)
}
if(colpar=="Black-White") {
   heatcol<-colorRampPalette(c("Black", "LightGrey"))(32)
}

# Plotting
par(mar=c(7, 4, 4, 2) + 0.1)
bitmap(file="heatmap.png", width=w/72, height=h/72)
heatmap(x=as.matrix(dat2), Rowv=clustg, Colv=clustc, col=heatcol, margins=c(10,5), labCol=gsub(" ", "", phenodata$description))
dev.off()
