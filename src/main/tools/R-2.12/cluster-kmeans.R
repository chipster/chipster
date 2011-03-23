# ANALYSIS Clustering/K-Means (K-means clustering of genes. Divides the genes in the selected data set into a specified 
# number of clusters.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT kmeans.tsv, kmeans.pdf
# PARAMETER number.of.clusters INTEGER FROM 2 TO 100000 DEFAULT 5 (Number of clusters)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the resampling image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the resampling image)


# K-means clustering
# JTT 22.6.2006

# Parameter settings (default) for testing purposes
#number.of.clusters<-10
#image.width<-600
#image.height<-600

# Renaming variables
k<-number.of.clusters
w<-image.width
h<-image.height

# Loads the data file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Sanity checks
if(k>=nrow(dat2)) {
   stop("You have specified as many clusters as there are genes! Use a smaller number of clusters.")
}

# Calculates the K-means clustering result
# Needs a parameter k, the number of clusters
# k<-c(5)
km<-kmeans(dat2, k, iter.max=100000, nstart=min(10, nrow(dat2)))

# Plotting the clustering
max.dat2<-max(dat2)
min.dat2<-min(dat2)
bitmap(file="kmeans.pdf", width=w/72, height=h/72)
par(mfrow=c(ceiling(sqrt(k)), ceiling(sqrt(k))))
for(i in 1:k) {
   matplot(t(dat2[km$cluster==i,]), type="l", main=paste("cluster:", i), ylab="log expression", col=1, lty=1, ylim=c(min.dat2, max.dat2))
}
dev.off()

# Writing a table
# Creates a table with column giving the cluster membership
dat2<-data.frame(dat, cluster=km$cluster)
write.table(dat2, "kmeans.tsv", sep="\t", row.names=T, col.names=T, quote=F)
