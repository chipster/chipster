# ANALYSIS Utilities/"Average replicate chips" (Calculates an arithmetic average of gene expression levels for replicate 
# chips.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT average-replicates.tsv
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to average)


# Two-group parametric and non-parametric tests
# JTT 30.7.2007

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Test needs a parameter "groups" that specifies the grouping of the samples
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
groups<-phenodata[,grep(column, colnames(phenodata))]

# Sanity checks
if(length(unique(groups))==1) {
   stop("You do not have any replicates to average!")
}

# Calculating averages
columnnames<-c()
dat3<-matrix(nrow=nrow(dat2), ncol=length(unique(groups)), NA)
for(i in 1:length(unique(groups))) {
   dat3[,i]<-rowSums(dat2[,which(groups==i)])
   columnnames<-c(columnnames, paste("group", i, sep=""))
}
columnnames<-paste("chip.", columnnames, sep="")
colnames(dat3)<-columnnames
rownames(dat3)<-rownames(dat2)

# Saving the results
write.table(data.frame(dat3), file="average-replicates.tsv", sep="\t", row.names=T, col.names=T, quote=F)
