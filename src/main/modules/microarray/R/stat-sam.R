# TOOL stat-sam.R: SAM (SAM analysis for one or more groups. You probably need to run this tool several times in order to get an idea of the results, since the results are reported for one delta-value at a time.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT sam.tsv: sam.tsv 
# OUTPUT sam.png: sam.png 
# OUTPUT sam-delta.pdf: sam-delta.pdf 
# PARAMETER column: "Column" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER p.value.threshold: "p-value threshold" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (FDR cut-off for significant results)
# PARAMETER random.number: "Random number" TYPE INTEGER FROM 1 TO 1000 DEFAULT 1 (Random number)
# PARAMETER number.of.delta: "Number of delta" TYPE INTEGER FROM 1 TO 100 DEFAULT 10 (How many different delta values are used)

# SAM analysis

# Parameter settings (default) for testing purposes
#column<-"group"
#p.value.adjustment.method<-"BH"
#p.value.threshold<-0.05
#random.number<-1
#number.of.delta<-10
#delta.to.plot<-1

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
samout<-sam(dat2, groups, rand=random.number, control=samControl(n.delta=number.of.delta))
delta_res <- findDelta(samout, fdr = p.value.threshold)
if(is.null(delta_res)) {
	write("No significant genes found *with this FDR cutoff*!", file="sam.tsv")

	bitmap(file="sam.png", width=600/72, height=600/72)
	plot(samout)
	dev.off()

	pdf(file="sam-delta.pdf", width=600/72, height=600/72)
	plot(samout, samout@mat.fdr[1,1])
	dev.off()
} else {
	if(length(delta_res) == 3) {
		delta_val <- delta_res[1]
	} else {
		delta_val <- delta_res[2,1]
	}

	# Plots the results
	bitmap(file="sam.png", width=600/72, height=600/72)
	plot(samout)
	dev.off()

	pdf(file="sam-delta.pdf", width=600/72, height=600/72)
	plot(samout, delta_val)
	dev.off()

	samsum<-summary(samout,delta_val)
	dat3<-dat2[samsum@row.sig.genes,]
	dat.sig<-samsum@mat.sig[,c(1,5,6)]
	dat.sig<-dat.sig[order(dat.sig$Row),]
	dat3<-data.frame(dat3, p.adjusted=dat.sig$q.value, FoldChange=dat.sig$R.fold)
	write.table(dat3, file="sam.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

