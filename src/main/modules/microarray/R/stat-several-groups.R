# TOOL stat-several-groups.R: "Several groups tests" (Tests for comparing the mean gene expression of several groups. Other than empiricalBayes might be slow, if run on unfiltered data.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT multiple-sample.tsv: multiple-sample.tsv 
# PARAMETER column: "Column" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER test: "Test" TYPE [empiricalBayes: "empirical Bayes", ANOVA: ANOVA, Kruskal-Wallis: Kruskal-Wallis] DEFAULT empiricalBayes (Test type)
# PARAMETER p.value.adjustment.method: "p-value adjustment method" TYPE [none: none, Bonferroni: Bonferroni, Holm: Holm, Hochberg: Hochberg, BH: BH, BY: BY] DEFAULT BH (Multiple testing correction method)
# PARAMETER p.value.threshold: "p-value threshold" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)

# column<-"group"
# test<-"empiricalBayes"
# p.value.adjustment.method <-"none"
# p.value.threshold <-0.05
# use.simple.analysis<-"no"

# JTT 4.7.2006 Created Several groups parametric and non-parametric tests
# MK 24.06.2013 
# EK 2.3.2014 Phenodata groups are always considered as factors

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
groups<-as.factor(phenodata[,pmatch(column,colnames(phenodata))])

# Sanity checks
if(ncol(dat2)<4) {
   stop("You need to have at least four chips to run this analysis")
}

# Empirical Bayes
if(meth=="empiricalBayes") {
  	library(limma)
  	design<-model.matrix(~groups)
   	fit<-lmFit(dat2, design)
   	fit<-eBayes(fit)

      tab <- toptable(fit, coef=2, number=nrow(fit), adjust.method=adj.method, sort.by="none")
      rows <- which(tab$adj.P.Val<=p.cut)
      p <- tab$adj.P.Val[rows]
      M <- tab$logFC[rows]
      dat <- data.frame(dat[rows,], p.adjusted=round(p, digits=6), FC=M);
      dat <- dat[order(dat$p.adjusted, decreasing=F),]

   	write.table(dat, file="multiple-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
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
	if (adj.method %in% c("Bonferroni", "Holm", "Hochberg")) {
		p.adjusted <- p.adjust(p.raw, method=tolower(adj.method))
	} else {
		p.adjusted <- p.adjust(p.raw, method=adj.method)		
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
	if (adj.method %in% c("Bonferroni", "Holm", "Hochberg")) {
		p.adjusted <- p.adjust(p.raw, method=tolower(adj.method))
	} else {
		p.adjusted <- p.adjust(p.raw, method=adj.method)		
	}

	dat<-dat[p.adjusted<=p.cut,]   
   	p.adjusted<-p.adjusted[p.adjusted<=p.cut]
   	write.table(data.frame(dat, p.adjusted=round(p.adjusted, digits=6)), file="multiple-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}