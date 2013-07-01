# TOOL cluster-hierarchical.R: Hierarchical (Hierarchical clustering of genes or chips. Specify the distance measure and the clustering method. The clustering result can be validated using bootstrapping. Validation is computationally very expensive, and works for approximately less than 100 genes. Clustering is done using the function hcluster in which the parameter correlation envokes computation of pearson type of distances.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT hc.tre: hc.tre 
# OUTPUT resample.pdf: resample.pdf 
# PARAMETER cluster: cluster TYPE [genes: genes, chips: chips] DEFAULT genes (What to cluster)
# PARAMETER distance.method: distance.method TYPE [euclidian: euclidian, manhattan: manhattan, pearson: pearson, spearman: spearman] DEFAULT pearson (Distance measure)
# PARAMETER tree.method: tree.method TYPE [single: single, average: average, complete: complete, ward: ward] DEFAULT average (Clustering method)
# PARAMETER resampling: resampling TYPE [none: none, bootstrap: bootstrap] DEFAULT none (Validation)
# PARAMETER number.of.replicates: number.of.replicates TYPE INTEGER FROM 1 TO 10000 DEFAULT 1000 (Number of pseudoreplicates to use in validation)
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the resampling image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the resampling image)

# Hierarchical clustering, permutation testing 
# JTT 22.6.2006
#
# MG 25.11.2010
# Increased the gene/sample limit to 20000

# Parameter settings (default) for testing purposes
#cluster<-c("chips")
#distance.method<-c("euclidian")
#tree.method<-c("ward")
#resampling<-c("none")
#number.of.replicates<-c(100)
#image.width<-600
#image.height<-600


# Renaming variables
margin<-cluster
distmeth<-distance.method
treemeth<-tree.method
doresample<-resampling
perms<-number.of.replicates
w<-image.width
h<-image.height
do.sample<-c("all")

if(distmeth=="pearson") {
   distmeth<-"correlation"
}

# Loads the libraries
library(amap)
library(ape)

# Loads the data file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Check that resampling is not applied to dataset larger than 1000 genes
if (resampling == "bootstrap" && nrow(dat)>1000) {
	stop("CHIPSTER-NOTE: Bootstrap resampling on datasets larger than 1000 genes is not possible due to computing time limitations. Please note that you can run hierarchical clustering on datasets including up to 20000 genes, provided the resampling option is turned off.")
}

# Takes a sample of the data
# PARAMETER do.sample [25, 50, 75, 100, 125, 150, 175, 200, all] DEFAULT 100 (If the data is large, should the image be generated from a specified number of random genes)
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

if (nrow(dat2) > 20000) {
  stop("Hierarchical clustering can be run on maximum 20000 of genes/samples");
}

# Tree calculation, no resampling
if(doresample=="none") {
   clust<-hcluster(x=dat2, method=distmeth, link=treemeth)
   phylo.clust<-as.phylo(clust)
   write.tree(phylo.clust, "hc.tre")
   # write.table(data.frame(dat, order=clust$order, mergex=c(clust$merge[,1], 0), mergey=c(clust$merge[,2], 0), height=c(clust$height, 0)), "hc.tre", sep="\t", row.names=T, col.names=T, quote=F)
}

# Tree calculation, with resampling
if(doresample=="none"){
   pdf(file="resample.pdf", width=w/72, height=h/72)
   plot(1, 1, col=0)
   text(1, 1, "This is a dummy image.", col=1)
   text(1, 0.9, "To generate a real image, turn on the bootstrapping option.", col=1)
   text(1, 0.8, "Note, that bootstrapping is not necessary, and might take a very long time to complete.", col=1)
   dev.off()
}
if(doresample=="bootstrap"){
   library(pvclust)
   clust<-hcluster(x=dat2, method=distmeth, link=treemeth)
   phylo.clust<-as.phylo(clust)
   write.tree(phylo.clust, "hc.tre")
   if(distmeth=="spearman" | distmeth=="manhattan") {
      stop("Resampling can only be applied if Pearson correlation or euclidian distance is used! If appropriate, change the settings accordingly.")
   }
   dat2<-t(dat2)
   pv.clust<-pvclust(dat2, method.dist=distmeth, method.hclust=treemeth, nboot=perms, r=1)
   pdf(file="resample.pdf", width=w/72, height=h/72)
   plot(pv.clust, cex.pv=0.75, font.pv=0.75, cex=0.75)
   dev.off()
}