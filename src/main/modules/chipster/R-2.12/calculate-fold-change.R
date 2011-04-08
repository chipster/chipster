# ANALYSIS Utilities/"Calculate fold change" (Calculates an arithmetic average of gene expression for replicate chips.
# Then calculates a ratio of the averages. Works only if you have exactly two groups of samples.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT fold-change.tsv
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to average.)


# Two-group parametric and non-parametric tests
# JTT 30.7.2007

# Parameter settings (default) for testing purposes
#column<-c("group")

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Test needs a parameter "groups" that specifies the grouping of the samples
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
groups<-phenodata[,pmatch(column,colnames(phenodata))]

# Sanity checks
if(length(unique(groups))==1) {
   stop("You do not have any replicates to average!")
}
if(length(unique(groups))>2) {
   stop("You have more than two groups! I don't know how to calculate fold change.")
}

# Calculating averages
columnnames<-c()
dat3<-matrix(nrow=nrow(dat2), ncol=length(unique(groups)), NA)
for(i in 1:length(unique(groups))) {
   dat3[,i]<-rowSums(data.frame(dat2[,which(groups==i)]))/ncol(data.frame(dat2[,which(groups==i)]))
   columnnames<-c(columnnames, paste("group", i, sep=""))
}
columnnames<-paste("chip.", columnnames, sep="")
colnames(dat3)<-columnnames
rownames(dat3)<-rownames(dat2)

# Calculating the fold change
# Treatment divided by the control
FD<-dat3[,2]-dat3[,1]

# Saving the results
write.table(data.frame(dat, logFC=FD), file="fold-change.tsv", sep="\t", row.names=T, col.names=T, quote=F)
