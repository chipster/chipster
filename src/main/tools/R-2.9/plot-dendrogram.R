# ANALYSIS Visualisation/"Dendrogram" (Creates a dendrogram of samples using normalized data with Pearson correlation and 
# average linkage method. The branches of the tree are colored according to the selected number of groups.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT dendrogram-color.png, dendrogram-bw.png
# PARAMETER cluster [genes, chips] DEFAULT chips (What to cluster)
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to color next to the tree)
# PARAMETER number.of.groups INTEGER FROM 2 TO 20 DEFAULT 2 (How many groups to color to the tree)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Dendrogram
# JTT 3.10.2007
#
# MG 25.11.2010
# Increased the gene/sample limit to 20000

# Renaming variables
w<-image.width
h<-image.height
gr<-number.of.groups 
margin<-cluster

# Loading the libraries
library(fpc)
library(A2R)
library(amap)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Loads the phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
groups<-phenodata[,pmatch(column,colnames(phenodata))]

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]
colnames(dat2)<-gsub(" ", "", phenodata$description)

# Manipulate data depending on what to cluster
if(margin=="chips") {
   dat2<-t(dat2)        
}

if (nrow(dat2) > 20000) {
  stop("Hierarchical clustering dendrogram can be plotted on maximum 2000 of genes/samples");
}

# Does the clustering
clust<-hcluster(x=dat2, method="pearson", link="average")
ct<-cutree(clust, gr)

# Plotting
bitmap(file="dendrogram-color.png", width=w/72, height=h/72)
if(margin=="chips") {
   A2Rplot(clust, k=gr, fact.sup=groups)        
} else {
   A2Rplot(clust, k=gr)        
}
dev.off()

bitmap(file="dendrogram-bw.png", width=w/72, height=h/72)
plot(clust)
dev.off()
