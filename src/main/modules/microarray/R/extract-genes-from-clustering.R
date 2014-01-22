# TOOL extract-genes-from-clustering.R: "Extract genes from clustering" (Extracts genes from a clustering result. Specify the cluster you want to extract the genes from.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT extract.tsv: extract.tsv 
# PARAMETER cluster.number: cluster.number TYPE INTEGER FROM 1 TO 1000 DEFAULT 1 (Cluster number to extract genes from.)


# Extracts genes from the clustering result for a specific cluster
# JTT 15.10.2007

# Loads the data file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Sanity checks
if(is.null(dat$cluster)==TRUE) {
   stop("You haven't clustered your data, so there are no genes to extract! Cluster the data first.")
}
if(cluster.number>max(dat$cluster)) {
   stop("The cluster number is larger than the number of cluster in the current data! Use a smaller cluster number.")
}

# Extracting the genes
dat2<-dat[which(dat$cluster==cluster.number),]

# Writing the data to disk
write.table(data.frame(dat2), file="extract.tsv", sep="\t", row.names=T, col.names=T, quote=F)

