# ANALYSIS Statistics/"SAM" (SAM analysis for one or more groups. You probably need to run this tool several
# times in order to get an idea of the results, since the results are reported for one delta-value at a time.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT sam.tsv, sam.png, sam-delta.png
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER p.value.adjustment.method [none, Bonferroni, BH] DEFAULT BH (Multiple testing correction method)
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER random.number INTEGER FROM 1 to 1000 DEFAULT 1 (Random number)
# PARAMETER number.of.delta INTEGER FROM 1 TO 100 DEFAULT 10 (How many different delta values are used)
# PARAMETER delta.to.plot DECIMAL FROM 0.1 TO 100 DEFAULT 1 (Which delta value to plot or to use for reporting)


# SAM analysis

#Loading the libraries
library(multtest)
library(siggenes)

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
if(column=="empty") {
   print("You haven't selected a phenodata column! Using group instead.")
   groups<-phenodata$group
}

# Runs the test
samout<-sam(dat2, groups, rand=random.number, n.delta=number.of.delta)

# Plots the results
bitmap(file="sam.png", width=600/72, height=600/72)
plot(samout)
dev.off()

bitmap(file="sam-delta.png", width=600/72, height=600/72)
plot(samout, delta.to.plot)
dev.off()

# Generates report
samsum<-summary(samout,delta.to.plot)
if(nrow(samsum@mat.sig)>=1) {
   dat3<-dat2[samsum@row.sig.genes,]
   dat.sig<-samsum@mat.sig[,c(1,5,6)]
   dat.sig<-dat.sig[order(dat.sig$Row),]
   dat3<-data.frame(dat3, p.adjusted=dat.sig$q.value, FoldChange=dat.sig$R.fold)
   write.table(dat3, file="sam.tsv", sep="\t", row.names=T, col.names=T, quote=F)
} else {
   write("No significant genes found *with this delta value*!", file="sam.tsv")
}

