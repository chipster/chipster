# Hierarchical clustering, two-way
# JTT, 14.7.2005

# Some things here are still broken, and should be fixed...
# Do we want to get out also the distance matrix or simply the tree?

# Loads the normalized data
# This needs a parameter "file" that tells which data file should be loaded.
# Datafile is produced by some statistical tool (currently, only stat-anova.R and filter-pvalue.R)
# file<-c("rma.txt")
# data<-read.table(file, header=T, sep="\t")
file<-c("sd-filter.txt", header=T, sep="\t")
data<-as.matrix(read.table(file, header=T, sep="\t"))

# Distance calculation
# This needs a parameter (distmeth) specifying the distance measure
# distmeth can be any of "euclidean", "maximum", "manhattan", "canberra", or "binary".
distmeth<-c("euclidian")
D<-dist(t(data), method=distmeth)

# Drawing the tree
# This needs a parameter (treemeth) specifying the tree drawing style
# treemeth can be any of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
treemeth<-c("average")
hc<-hclust(D, method=treemeth)

# Drawing a heatmap
cwd=getwd()
jpeg(file.path(cwd, "heatmap.jpg"), quality=100)
heatmap(data, distfun = dist, hclustfun = hclust)
dev.off()