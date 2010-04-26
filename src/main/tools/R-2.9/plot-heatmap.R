# ANALYSIS Visualisation/"Heatmap" (Draws a heatmap using Pearson correlation and average linkage.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT heatmap.png
# PARAMETER coloring.scheme [Green-Red, Blue-Yellow, Black-White] DEFAULT Blue-Yellow (Coloring scheme for the SOM map)
# PARAMETER cluster.samples.only [yes, no] DEFAULT no (Disables clustering on the genes. This option is convenient if you
# want to retain a predefined gene order or make a sample clustering heatmap with more than 10000 genes)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Heatmap
# JTT 3.10.2007
#
# modified by MG, 21.4.2010, to increase the gene/sample limit and adjust column width to number of genes
# and have the option to cluster samples only

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

# Sanity checks
if (nrow(dat2) > 10000 & cluster.samples.only=="no") {
	stop("CHIPSTER-NOTE: Hierarchical clustering can be run on a maximum of 10000 genes");
}
if (ncol(dat2) > 10000) {
	stop("CHIPSTER-NOTE: Hierarchical clustering can be run on a maximum of 10000 samples");
}

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
colnames(dat2)<-gsub(" ", "", phenodata$description)

# Does the clustering
if (cluster.samples.only=="no") {
	clustg<-as.dendrogram(hcluster(x=dat2, method="pearson", link="average"))
}
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

	# set up column margin according to number of genes
	if (nrow(dat2)>200) {
		column_margin <- 5
	}
	if (nrow(dat2)<=200 & nrow(dat2)>50) {
		column_margin <- 10
	}
	if (nrow(dat2)<=50) {
		column_margin <- 15
	}
	
bitmap(file="heatmap.png", width=w/72, height=h/72)
if (cluster.samples.only=="no") {
	heatmap(x=as.matrix(dat2), Rowv=clustg, Colv=clustc, col=heatcol, margins=c(15, column_margin), labCol=gsub(" ", "", phenodata$description))
}
if (cluster.samples.only=="yes") {
	heatmap(x=as.matrix(dat2), Rowv=NA, Colv=clustc, col=heatcol, margins=c(15, column_margin), labCol=gsub(" ", "", phenodata$description))
}
dev.off()
