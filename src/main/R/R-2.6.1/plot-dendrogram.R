# ANALYSIS Visualisation/"Dendrogram" (Creates a dendrogram of samples using normalized data with Pearson correlation and 
# average linkage method. The branches of the tree are colored according to the selected number of groups.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT dendrogram-color.png, dendrogram-bw.png
# PARAMETER cluster [genes, chips] DEFAULT chips (What to cluster)
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to color next to the tree)
# PARAMETER number.of.groups INTEGER FROM 1 TO 20 DEFAULT 2 (How many groups to color to the tree)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Dendrogram
# JTT 3.10.2007

# Renaming variables
w<-image.width
h<-image.height
gr<-number.of.groups 
margin<-cluster
do.sample<-c("all")

# Loading the libraries
library(fpc)
library(A2R)
library(amap)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Loads the phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
groups<-phenodata[,grep(column, colnames(phenodata))]

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

# Manipulate data depending on what to cluster
if(margin=="chips") {
   dat2<-t(dat2)        
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