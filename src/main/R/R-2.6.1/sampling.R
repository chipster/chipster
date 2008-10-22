# ANALYSIS Utilities/"Random sampling" (Generates a random sample of the data.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT sample.tsv
# PARAMETER what.to.sample [genes, chips] DEFAULT genes (Sample genes or chips)
# PARAMETER sample.size DECIMAL FROM 0 TO 100 DEFAULT 10 (Percentage of the data to sample)


# Random sampling
# JTT 20.10.2008

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Sanity checks
if(ncol(dat)<2 | nrow(dat)<2) {
   stop("You need to have at least 2 data points from which to sample!")
}

# Sampling
if(what.to.sample=="genes") {
   v<-sample(1:nrow(dat), floor(nrow(dat)*sample.size/100))
   dat2<-dat[v,]
}

if(what.to.sample=="chips") {
   v<-sample(1:ncol(dat), floor(ncol(dat)*sample.size/100))
   dat2<-dat[,v]
}

# Writing data to disk
write.table(data.frame(dat2), file="sample.tsv", col.names=T, quote=F, sep="\t", row.names=T)
