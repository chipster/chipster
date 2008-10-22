# ANALYSIS "Clustering"/"K-Means" (K-means clustering, clusters genes)
# INPUT GENE_EXPRS normalised.txt OUTPUT kmeans.txt
# PARAMETER k INTEGER FROM 0 TO 100 DEFAULT 5 (number of clusters)

# JTT, 14.7.2005

# Loads the normalized data
# This needs a parameter "file" that tells which data file should be loaded.
# Datafile is produced by some statistical tool (currently, only stat-anova.R)
# At the moment datafile can be only be "p-anova.txt"
file<-c("normalised.txt")
data<-read.table(file, header=T, sep="\t")

# Calculates the K-means clustering result
# This needs a parameter k, the number of clusters (inputted by the user)
# In the result, $cluster gives the cluster the gene belongs to,
# $withinss gives the similarity of the genes, and $size the number of
# genes in each cluster
#k<-c(5)
km<-kmeans(data, k, iter.max=1000)

# Modifying the data for output
# Now the data has one more column giving the class membership (cl)
cl<-km$cluster
data.kmeans<-data.frame(data,cl)

# Results to be shown for the user:
# Clustering result (drawn from "kmeans.txt"), sizes of the clusters,
# and withinness of the clusters (clusterwise). Withinness is calculated 
# here. 
clusterwise.withinness<-km$withinss/km$size

# Saving the results
write.table(data.kmeans, "kmeans.txt", sep="\t", row.names=T, col.names=F) # (akallio) changed col.names=F

# This would produce a nice image of the first two chips against each other
# and genes colored by the K-means cluster membership
# plot(data[,1], data[,2], col=K$cluster)
# points(km$centers, col = 1:2, pch = 8)

# This would save the crude K-means results
# sink("kmeans.txt")
# K
# sink()