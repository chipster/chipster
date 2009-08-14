# ANALYSIS Statistics/"Several groups tests" (Tests for comparing the mean gene expression of several groups.
# Other than empiricalBayes might be slow, if run on unfiltered data.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT multiple-sample.tsv
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER test [empiricalBayes, ANOVA, Kruskal-Wallis] DEFAULT empiricalBayes (Test type)
# PARAMETER p.value.adjustment.method [none, Bonferroni, Holm, Hochberg, BH, BY] DEFAULT BH (Multiple testing correction method)
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER use.simple.analysis [yes, no] DEFAULT no (Use simple analysis: if no significant genes are found, a hundred most significant genes are reported)


# column<-"group"
# test<-"empiricalBayes"
# p.value.adjustment.method <-"none"
# p.value.threshold <-0.05
# use.simple.analysis<-"no"


# Several groups parametric and non-parametric tests
# JTT 4.7.2006

# Loads the libraries
library(multtest)

# Renaming variables
meth<-test
if(test=="empiricalBayes" & (p.value.adjustment.method!="BH" & p.value.adjustment.method!="BY") ) {
   adj.method<-tolower(p.value.adjustment.method)
} else {
   adj.method<-p.value.adjustment.method
}
p.cut<-p.value.threshold

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
if(ncol(dat2)<4) {
   stop("You need to have at least four chips to run this analysis")
}

# Testing

if(use.simple.analysis=="no") {
# Empirical Bayes
if(meth=="empiricalBayes") {
   library(limma)
   design<-model.matrix(~groups)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
   tab<-toptable(fit, coef=2, number=nrow(fit), adjust.method=adj.method)
   rows<-as.numeric(row.names(tab))
   rows<-rows[tab$adj.P.Val<=p.cut]
   p<-tab$adj.P.Val[tab$adj.P.Val<=p.cut]
   dat<-dat[rows,]
   write.table(data.frame(dat, p.adjusted=round(p, digits=6)), file="multiple-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# ANOVA test
if(meth=="ANOVA") {
   p<-c()
   len<-length(dat2[,1])
   for(i in 1:len) {
      p<-c(p, na.omit((anova(lm(t(dat2)[,i]~groups)))$Pr)[1])
   }
   p.raw<-p
   # Multiple testing correction
   if(adj.method=="none") {
      p.adjusted<-p.raw
   }
   if(adj.method=="Bonferroni" | adj.method=="BH") {
      p.adjusted<-mt.rawp2adjp(p.raw, adj.method)
      p.adjusted<-p.adjusted$adjp[order(p.adjusted$index),][,2]
   }
   dat<-dat[p.adjusted<=p.cut,]   
   p.adjusted<-p.adjusted[p.adjusted<=p.cut]
   write.table(data.frame(dat, p.adjusted=round(p.adjusted, digits=6)), file="multiple-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# Kruskal-Wallis test
if(meth=="Kruskal-Wallis") {
   p<-c()
   for(i in 1:nrow(dat2)) {
      p<-c(p, kruskal.test(dat2[i,], groups))
   }
   p.raw<-p
   if(adj.method=="none") {
      p.adjusted<-p.raw
   }
   if(adj.method=="Bonferroni" | adj.method=="BH") {
      p.adjusted<-mt.rawp2adjp(p.raw, adj.method)
      p.adjusted<-p.adjusted$adjp[order(p.adjusted$index),][,2]
   }
   dat<-dat[p.adjusted<=p.cut,]   
   p.adjusted<-p.adjusted[p.adjusted<=p.cut]
   write.table(data.frame(dat, p.adjusted=round(p.adjusted, digits=6)), file="multiple-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}
}

if(use.simple.analysis=="yes") {
# Empirical Bayes
if(meth=="empiricalBayes") {
   library(limma)
   design<-model.matrix(~groups)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
   if(adj.method=="none") {
      tab<-toptable(fit, coef=2, number=nrow(fit), adjust.method="none")
      p.adjusted<-tab$adj.P.Val
      p<-tab$P.Value
   }
   if(adj.method=="BH") {
      tab<-toptable(fit, coef=2, number=nrow(fit), adjust.method="BH")
      p.adjusted<-tab$adj.P.Val
      p<-tab$P.Value
   }
   p.adjusted[p.adjusted>p.cut]<-NA
}

# ANOVA test
if(meth=="ANOVA") {
   p<-c()
   len<-length(dat2[,1])
   for(i in 1:len) {
      p<-c(p, na.omit((anova(lm(t(dat2)[,i]~groups)))$Pr)[1])
   }
   p.raw<-p
   if(adj.method=="none") {
      p.adjusted<-p.raw
   }
   if(adj.method=="Bonferroni" | adj.method=="BH") {
      p.adjusted<-mt.rawp2adjp(p.raw, adj.method)
      p.adjusted<-p.adjusted$adjp[order(p.adjusted$index),][,2]
   }
   p.adjusted[p.adjusted>p.cut]<-NA
}

# Kruskal-Wallis test
if(meth=="Kruskal-Wallis") {
   p<-c()
   for(i in 1:nrow(dat2)) {
      p<-c(p, kruskal.test(dat2[i,], groups))
   }
   p.raw<-p
   if(adj.method=="none") {
      p.adjusted<-p.raw
   }
   if(adj.method=="Bonferroni" | adj.method=="BH") {
      p.adjusted<-mt.rawp2adjp(p.raw, adj.method)
      p.adjusted<-p.adjusted$adjp[order(p.adjusted$index),][,2]
   }
   p.adjusted[p.adjusted>p.cut]<-NA
}

# Writing the data to disk
if(nrow(na.omit(data.frame(dat, p.adjusted=round(p.adjusted, digits=6))))==0) {
   dat<-data.frame(dat[order(p),], p.raw=round(p[order(p)], digits=6))
   dat<-dat[1:100,]
} else {
   dat<-na.omit(data.frame(dat, p.adjusted=round(p.adjusted, digits=6)))
}
write.table(dat, file="multiple-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}
