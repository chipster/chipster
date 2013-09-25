# TOOL plot-heatmap.R: Heatmap (Draws a heatmap using Pearson correlation and average linkage. Clustering is done using the function hcluster in which the parameter correlation envokes computation of pearson type of distances.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT heatmap.pdf: heatmap.pdf 
# PARAMETER coloring.scheme: coloring.scheme TYPE [Green-Red: Green-Red, Blue-Yellow: Blue-Yellow, Black-White: Black-White] DEFAULT Blue-Yellow (Coloring scheme for the SOM map)
# PARAMETER cluster.samples.only: cluster.samples.only TYPE [yes: yes, no: no] DEFAULT no (Disables clustering on the genes. This option is convenient if you want to retain a predefined gene order or make a sample clustering heatmap with more than 10000 genes)
# PARAMETER hm.scale: hm.scale TYPE [default: default, none: none, row: row, column: column] DEFAULT default (Indicating if the values should be centered and scaled in either the row direction or the column direction, or none.)
# PARAMETER distance: distance TYPE [euclidean: euclidean, maximum: maximum, manhattan: manhattan, canberra: canberra, binary: binary, pearson: pearson, spearman: spearman, kendall: kendall] DEFAULT pearson (The correlation distance measure to be used for clustering.)
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Heatmap
# JTT 3.10.2007
#
# modified by MG, 21.4.2010, to increase the gene/sample limit and adjust column width to number of genes
# and have the option to cluster samples only
#
# MG 25.11.2010
# Increased the gene/sample limit to 20000
#
# OH 10.10.2012
# removed dependency library("amap")
# added library("MKmisc)
# added parameter hm.scale, see parameter scale of heatmap function, R: help(heatmap)
# added parameter distance, see parameter distfun of heatmap function, to allow for different distance measurements
#
# OH 11.06.2013
# switching back from function hclust (used by function heatmap) to pre-calc dendrogram with function hcluster:
#   visualization "plot-dendrogram.R" produced a different dendrogram than plot-heatmap.R with same parameters
#   because of usage of different functions hcluster and hclust.
#


# Parameter settings (default) for testing purposes
#coloring.scheme<-"Green-Red"
#cluster.samples.only<-"no"
#image.width<-600
#image.height<-600

# Renaming variables
colpar<-coloring.scheme
w<-image.width
h<-image.height

s<-hm.scale
d<-distance

if(d=="pearson") {
   d<-"correlation"
}

# Loading the libraries
library(amap)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Sanity checks
if (nrow(dat2) > 20000 & cluster.samples.only=="no") {
	stop("CHIPSTER-NOTE: Hierarchical clustering can be run on a maximum of 20000 genes");
}
if (ncol(dat2) > 20000) {
	stop("CHIPSTER-NOTE: Hierarchical clustering can be run on a maximum of 20000 samples");
}

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
colnames(dat2)<-gsub(" ", "", phenodata$description)

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

# Does the clustering
#if (cluster.samples.only=="no") {
#	clustg<-as.dendrogram(hcluster(x=dat2, method="pearson", link="average"))
#}
#clustc<-as.dendrogram(hcluster(x=t(dat2), method="pearson", link="average"))

#library("MKmisc")
#if( d=="euclidean" || d=="maximum" || d=="manhattan" || d=="canberra" || d=="binary" || d=="minkowski" ) {
#	myCorDist=function(x){dist(x,method=d)}
#}
#if( d=="pearson" || d=="spearman" || d=="kendall" || d=="cosine" || d=="mcd" || d=="ogk" ) {
#	myCorDist=function(x){corDist(x,method=d)}
#}

if (cluster.samples.only=="no") {
	clustg<-as.dendrogram(hcluster(x=dat2, method=d, link="average"))
}
clustc<-as.dendrogram(hcluster(x=t(dat2), method=d, link="average"))

pdf(file="heatmap.pdf", width=w/72, height=h/72)
if (cluster.samples.only=="no") {
	if( hm.scale=="default" ) {
		heatmap(x=as.matrix(dat2), Rowv=clustg, Colv=clustc, col=heatcol, margins=c(15, column_margin), labCol=gsub(" ", "", phenodata$description))
		#heatmap(x=as.matrix(dat2), distfun=myCorDist, hclustfun=hclust, col=heatcol, margins=c(15, column_margin), labCol=gsub(" ", "", phenodata$description))
	} else {
		heatmap(x=as.matrix(dat2), Rowv=clustg, Colv=clustc, col=heatcol, scale=s, margins=c(15, column_margin), labCol=gsub(" ", "", phenodata$description))
		#heatmap(x=as.matrix(dat2), distfun=myCorDist, hclustfun=hclust, col=heatcol, scale=s, margins=c(15, column_margin), labCol=gsub(" ", "", phenodata$description))
	}
}
if (cluster.samples.only=="yes") {
	if( hm.scale=="default" ) {
		heatmap(x=as.matrix(dat2), Rowv=NA, Colv=clustc, col=heatcol, margins=c(15, column_margin), labCol=gsub(" ", "", phenodata$description))
		#heatmap(x=as.matrix(dat2), Rowv=NA, distfun=myCorDist, hclustfun=hclust, col=heatcol, margins=c(15, column_margin), labCol=gsub(" ", "", phenodata$description))
	} else {
		heatmap(x=as.matrix(dat2), Rowv=NA, Colv=clustc, col=heatcol, scale=s, margins=c(15, column_margin), labCol=gsub(" ", "", phenodata$description))
		#heatmap(x=as.matrix(dat2), Rowv=NA, distfun=myCorDist, hclustfun=hclust, col=heatcol, scale=s, margins=c(15, column_margin), labCol=gsub(" ", "", phenodata$description))
	}
}
dev.off()
