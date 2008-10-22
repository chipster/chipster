# ANALYSIS Visualisation/"Heatmap" (Draws a heatmap using Pearson correlation and average linkage.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT heatmap.png
# PARAMETER coloring.scheme [Green-Red, Blue-Yellow, Black-White] DEFAULT Blue-Yellow (Coloring scheme for the SOM map)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Heatmap
# JTT 3.10.2007

# Renaming variables
colpar<-coloring.scheme
w<-image.width
h<-image.height
do.sample<-c("all")

# Loading the libraries
library(amap)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]
colnames(dat2)<-gsub("chip.(.+)", "\\1", colnames(dat2)[grep("chip", colnames(dat2))])

# Takes a sample of the data
# PARAMETER do.sample [25, 50, 75, 100, 125, 150, 175, 200, all] DEFAULT all (If the data is large, should the image be generated from a specified number of random genes)
if(do.sample!="all") {
   sample.size<-as.numeric(do.sample)
   if(nrow(dat2)<sample.size) {
      do.sample<-"all"
   }
   if(do.sample=="25") {
      dat2<-dat2[sample(nrow(dat2), 25),]
   }
   if(do.sample=="50") {
      dat2<-dat2[sample(nrow(dat2), 50),]
   }
   if(do.sample=="75") {
      dat2<-dat2[sample(nrow(dat2), 75),]
   }
   if(do.sample=="100") {
      dat2<-dat2[sample(nrow(dat2), 100),]
   }
   if(do.sample=="125") {
      dat2<-dat2[sample(nrow(dat2), 125),]
   }
   if(do.sample=="150") {
      dat2<-dat2[sample(nrow(dat2), 150),]
   }
   if(do.sample=="175") {
      dat2<-dat2[sample(nrow(dat2), 175),]
   }
   if(do.sample=="200") {
      dat2<-dat2[sample(nrow(dat2), 200),]
   }
}

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
heatmap(x=as.matrix(dat2), Rowv=clustg, Colv=clustc, col=heatcol)
dev.off()
