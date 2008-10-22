# ANALYSIS Pathways/"*Boolean network" (Searches for the best Boolean network for the specified timeseries. You need to
# have a time column in your phenodata, otherwise the analysis will fail miserably.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT best-network.tsv
# PARAMETER number.of.inputs.per.fuction INTEGER FROM 1 TO 100 DEFAULT 2 (Number of inputs per function in the network)
# PARAMETER number.of.genes INTEGER FROM 1 TO 100 DEFAULT 50 (Number of genes used for network inference)
# PARAMETER discretization.method [median, cutoff] DEFAULT median (Discretization method)


# bestim
# JTT 13.11.2006

# Loading libraries
library(e1071)
library(GeneTS)
library(genefilter)

# Reading data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Loading sample data
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Renaming variables
k<-number.of.inputs.per.fuction
genes<-number.of.genes
disc.meth<-discretization.method

# How many replicates there are per time point?
if(length(phenodata$time)>length(levels(as.factor(phenodata$time)))) {
   repl<-c()
   for(i in 1:length(levels(as.factor(phenodata$time)))) {
      repl<-c(repl, sum(as.numeric(grep(phenodata$time[1], phenodata$time, value=T)==phenodata$time[1])))
   }
} else {
   repl<-rep(1, length(phenodata$time))
}

# Making a longitudinal object
dat3<-as.longitudinal(t(dat2), repeats=repl, time=phenodata$time)

# Replacing missing values
dat3<-t(na.omit(t(impute(dat3))))

# Scaling the data to the same mean
dat3<-t(genescale(t(dat3)))

# Finding periodically expressed genes
f<-fisher.g.test(dat3)
e<-fdr.estimate.eta0(f, method="conservative", diagnostic.plot=F)
p.adj<-fdr.control(f, Q = 0.05, eta0=e)$qvalues

# Getting the most significant genes (no. set by the user)
dat4<-dat3[,order(p.adj)[1:genes]]

# Recoding the expression with binary variables
if(disc.meth=="median") {
   dat5<-t(ifelse(dat4<=0, 0, 1))
}
if(disc.meth="cutoff") {
   dat5<-t(ifelse(dat4<=-1 | dat4 >=1, 1, 0))
}

# Recoding into state changes
# This part is tricky, and currently unfinished
#manh<-diag(as.matrix(dist(t(dat5), method="manhattan", upper=T))[,-1])
#manh<-ifelse(manh<median(manh), 0, 1)
#zeros<-grep(0, manh)
#ones<-grep(1, manh)
#for(i in 1:(length(zeros)-1)) {
#   if(zeros[i+1]>zeros[i]+1) {
#      print(i)
#   }
#}
#for(i in 1:(length(ones)-1)) {
#   if(ones[i+1]>ones[i]+1) {
#      print(i)
#   }
#}

# Creates input matrices and weights for bestim
X<-dat5[,-ncol(dat5)]
Y<-dat5[,-1]
W<-matrix(rep(1, ncol(X)), nrow=1, ncol=ncol(X))
unlink("X.tsv")
write(file="X.tsv", paste(nrow(X), "\t", ncol(X), sep=""), append=F)
write.table(x=X, file="X.tsv", sep="\t", row.names=F, col.names=F, quote=F, append=T)
unlink("Y.tsv")
write(file="Y.tsv", paste(nrow(Y), "\t", ncol(Y), sep=""), append=F)
write.table(x=Y, file="Y.tsv", sep="\t", row.names=F, col.names=F, quote=F, append=T)
unlink("W.tsv")
write(file="W.tsv", paste(nrow(W), "\t", ncol(W), sep=""), append=F)
write.table(x=W, file="W.tsv", sep="\t", row.names=F, col.names=F, quote=F, append=T)

# Runs the analysis
system(eval(paste("/v/solaris9/appl/molbio/sysbio/bestfit/bin/bnBestFit --X X.tsv --Y Y.tsv --w W.tsv --Xt X.tsv --Yt Y.tsv --Fhat fhat.tsv --Ehat ehat.tsv --Et et.tsv", " -k ", k, sep="")))

# Picking the best network
ehat<-read.table("ehat.tsv", sep="\t", header=F, skip=1)[,1:genes]
fhat<-read.table("fhat.tsv", sep="\t", header=F, skip=1)
fhat<-fhat[,-ncol(fhat)]
ehat.sums<-as.numeric(rowSums(ehat))
best<-fhat[which(ehat.sums==min(ehat.sums)):(which(ehat.sums==min(ehat.sums))+2^k-1),]
write.table(best, file="best-network.tsv", sep="\t", col.names=F, row.names=F, quote=F)
