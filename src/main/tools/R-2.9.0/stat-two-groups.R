# ANALYSIS Statistics/"Two groups tests" (Tests for comparing the mean gene expression of two groups. 
# LPE only works, if the whole data is used, i.e., the data should not be pre-filtered, if LPE is used. 
# Other than empiricalBayes might be slow, if run on unfiltered data.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT two-sample.tsv
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER test [empiricalBayes, fast-t-test, t-test, F-test, Mann-Whitney, LPE] DEFAULT empiricalBayes (Test type)
# PARAMETER p.value.adjustment.method [none, Bonferroni, Holm, Hochberg, BH, BY] DEFAULT BH (Multiple testing correction method)
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)


# Two-group parametric and non-parametric tests
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
if(length(unique(groups))==1 | length(unique(groups))>=3) {
   stop("You need to have exactly two groups to run this analysis")
}

# Testing

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
   M<-tab$logFC[tab$adj.P.Val<=p.cut]
   dat<-dat[rows,]
   write.table(data.frame(dat, p.adjusted=round(p, digits=6), FC=M), file="two-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# Fast T-test
if(meth=="fast-t-test") {
   fit1<-lm(t(dat2)~groups)
   p<-rep(NA, nrow(dat2))
   for(i in 1:nrow(dat2)) {
      sum(fit1$residuals[,i]^2)->sse
      sum((dat2[i,]-mean(as.numeric(dat2[i,])))^2)->sst
      r2<-1-(sse/sst)
      f<-r2/((1-r2)/(ncol(dat2)-1))
      p[i]<-1-pf(f, 1, (ncol(dat2)-1))
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
   write.table(data.frame(dat, p.adjusted=round(p.adjusted, digits=6)), file="two-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# T-test
if(meth=="t-test") {
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
   dat<-dat[p.adjusted<=p.cut,]   
   p.adjusted<-p.adjusted[p.adjusted<=p.cut]
   write.table(data.frame(dat, p.adjusted=round(p.adjusted, digits=6)), file="two-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# F-test
if(meth=="F-test") {
   p<-c()
   len<-nrow(dat2)
   dat2.1<-dat2[,groups==unique(groups)[1]]
   dat2.2<-dat2[,groups==unique(groups)[2]]
   for(i in 1:len) {
      p<-c(p, var.test(x=as.numeric(dat2.1[i,]), y=as.numeric(dat2.2[i,]))$p.value)
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
   write.table(data.frame(dat, p.adjusted=round(p.adjusted, digits=6)), file="two-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}


# Mann-Whitney test
if(meth=="Mann-Whitney") {
   dat3<-split(as.data.frame(t(dat2)), groups)
   # Split creates a list with two objects
   g1<-dat3$'1'
   g2<-dat3$'2'
   p<-c()
   for(i in 1:nrow(dat2)) {
      p<-c(p, wilcox.test(g1[,i], g2[,i], mu=0)$p.value)
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
   write.table(data.frame(dat, p.adjusted=round(p.adjusted, digits=6)), file="two-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

if(meth=="LPE") {
   library(LPE)
   group1<-dat2[,which(groups==1)]
   group2<-dat2[,which(groups==2)]
   g1.x<-baseOlig.error(group1)
   g2.x<-baseOlig.error(group2)
   lp<-data.frame(lpe(group1, group2, g1.x, g2.x, probe.set.name=row.names(dat2)))
   if(adj.method=="none") {
      x.location <- grep("^x", names(lp))
      y.location <- grep("^y", names(lp))
      x <- lp[, x.location]
      y <- lp[, y.location]
      pnorm.diff <- pnorm(lp$median.diff, mean = 0, sd = lp$pooled.std.dev)
      fdr <- 2 * apply(cbind(pnorm.diff, 1 - pnorm.diff), 1, min)
      fdr<-round(fdr, digits=4)
      dat2<-data.frame(dat, p.adjusted=fdr)
   }
   if(adj.method=="Bonferroni" | adj.method=="BH") {
      fdr<-fdr.adjust(lp, adjp=adj.method)
      dat2<-merge(dat,as.data.frame(round(fdr, digits=4)), by.x="row.names", by.y="row.names")
      dat2<-dat2[,-ncol(dat2)] # Removes the last columns that holds Z-test values
      names(dat2)[which(names(dat2)=="FDR")]<-"p.adjusted" # Renames "FDR" with "p.adjusted"
   }
   dat<-dat[p.adjusted<=p.cut,]   
   p.adjusted<-p.adjusted[p.adjusted<=p.cut]
   write.table(data.frame(dat, p.adjusted=round(p.adjusted, digits=6)), file="two-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}
