# ANALYSIS Statistics/"One sample tests" (Tests whether genes are differentially expressed, i.e., whether their expression
# is significantly different from the expected.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT one-sample.tsv
# PARAMETER test [t-test, Wilcoxon] DEFAULT t-test (Test type)
# PARAMETER scale.to.same.mean [yes, no] DEFAULT yes (Scale the data to the same mean before filtering)
# PARAMETER p.value.adjustment.method [none, Bonferroni, Holm, Hochberg, BH, BY] DEFAULT BH (Multiple testing correction method)
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER assumed.mean DECIMAL FROM -100000 TO 100000 DEFAULT 0 (Assumed mean of the data)


# One-sample parametric and non-parametric tests
# JTT 4.7.2006

# Load the libraries
library(multtest)
library(genefilter)

# Renaming variables
meth<-test
adj.method<-p.value.adjustment.method
p.cut<-p.value.threshold
mu<-assumed.mean

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Sanity checks
if(ncol(dat2)==1) {
   stop("You need to have at least two chips to run this analysis!")
}

# Scaling the data to the same mean
if(scale.to.same.mean=="yes") {
   scaled.dat2<-genescale(dat2)
} else {
   scaled.dat2<-dat2
}

# Testing
if(meth=="t-test") {
   library(genefilter)
   n<-ncol(scaled.dat2)
   d<-sqrt(n)
   s<-rowSds(scaled.dat2)/d
   m<-rowSums(scaled.dat2)/n
   t<-(m-mu)/s
   p.raw<-(1-pt(abs(t), df=n-1))*2
}

if(meth=="Wilcoxon") {
   p<-c()
   len<-length(scaled.dat2[,1])
   for(i in 1:len) {
      # Calculates the Wilcoxon test statistic for every row (gene)
      p<-c(p, wilcox.test(x=as.vector(as.numeric(scaled.dat2[i,])), mu=mu, alternative=c("two.sided"))[[3]])
   }
   p.raw<-p
}

# Multiple testing correction
if(adj.method=="none") {
   p.adjusted<-p.raw
}
if(adj.method!="none") {
   p.adjusted<-mt.rawp2adjp(p.raw, adj.method)
   p.adjusted<-p.adjusted$adjp[order(p.adjusted$index),][,2]
}

# Discarding p-values that do not fulfill the cutoff (p.cut)
dat<-dat[p.adjusted<=p.cut,]   
p.adjusted<-p.adjusted[p.adjusted<=p.cut]

# Writing out a table
write.table(data.frame(dat, p.adjusted=round(p.adjusted, digits=6)), file="one-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
