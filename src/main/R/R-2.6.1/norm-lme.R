# ANALYSIS Normalisation/"Random effects" (Removes the possible effect of random effects from the normalized expression
# values of every gene. This can be used for removing batch effects due to day, technician, etc. from the data. 
# You should have a column in your phenodata that defines the random effects groups, currently
# only two groups are allowed.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT normalized-lme.tsv
# PARAMETER column.groups METACOLUMN_SEL DEFAULT group (Phenodata column containing group effects)
# PARAMETER column.random METACOLUMN_SEL DEFAULT random (Phenodata column containing random effects groups)

# Linear Mixed Model
# 
# JTT 12.7.2006
# Modified to use nlme library on 19.10.2006

# Loads the libraries
library(nlme)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Test needs a parameter "groups" that specifies the grouping of the samples
# and a parameter random that specifies the random effect such as day or technician
groups<-phenodata[,grep(column.groups, colnames(phenodata))]
random<-phenodata[,grep(column.random, colnames(phenodata))]

# Sanity check
if(nrow(random)==0) {
   stop("You haven't specified any random effects! Please modify the phenodata before running this tool again.")
}
if(length(unique(random))==1) {
   stop("You only have one level in variable random. Specify at least two levels for the random effect.")
}

# Fits a linear mixed model for every gene, and saves residuals into a table
# Assumes no interaction between random and groups
dat3<-dat2
for(i in 1:nrow(dat3)) {
   gene<-as.vector(as.numeric(dat3[i,]))
   dat3[i,]<-resid(lme(fixed=gene~groups, random=gene~1|random, control=list(maxIter=10000), method="REML"))
}

# Writes a table of results
write.table(data.frame(round(dat3, digits=2), calls), file="normalized-lme.tsv", sep="\t", row.names=T, col.names=T, quote=F)
