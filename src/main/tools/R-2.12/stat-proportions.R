# Test for differences in proportions
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER test [Fisher, Chisquare, hypergeometric] DEFAULT Chisquare (Test type)
# PARAMETER p.value.adjustment.method [none, Bonferroni, Holm, Hochberg, BH, BY] DEFAULT BH (Multiple testing correction method)
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)


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
if(length(unique(groups))==1 | length(unique(groups))>=3) {
	stop("You need to have exactly two groups to run this analysis")
}

# Chi-square test
if(meth=="Chisquare") {
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



