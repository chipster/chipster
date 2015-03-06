# TOOL stat-hyperG-safe.R: "SAFE test for KEGG pathway enrichment" (Finds pathways from the Kyoto Encyclopedia of Genes and Genomes, KEGG, that are over- or under-represented in the selected gene list.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT safeplot.png: safeplot.png 
# OUTPUT safe.tsv: safe.tsv 
# PARAMETER column: Column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER p.value.threshold: "p-value threshold" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)
# PARAMETER minimum.category.size: "Minimum category size" TYPE INTEGER FROM 1 TO 100 DEFAULT 10 (Minimum size for categories to be evaluated)
# PARAMETER num.or.cat: "Phenodata type" TYPE [factor: factor, continuous: continuous] DEFAULT factor (Type of the phenodata column. Can either be factor or continuous)

# MG 5.11.2009: Accounting for changes in SparseM package regarding matrix.csr class
# MK 29.08.2013: Bug preventing ploting of plots without any genes in the shaded region fixed. Added ability to analyse factorial categories and numeric categories

# PARAMETER which.ontology [KEGG, GO] DEFAULT KEGG (Which ontology to use in the test?)
# Parameter settings (default) for testing purposes
#column<-"group"
#which.ontology<-"KEGG"
#p.value.threshold<-1
#image.width<-600
#image.height<-600
#minimum.category.size<-10

# Translating the variables
w<-image.width
h<-image.height

# Loads the libraries
library(SparseM)
library(safe)
library(multtest)

# Reads the chiptype from phenodata table
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
if(phenodata$chiptype[1]!="cDNA" | phenodata$chiptype[1]!="Illumina") {
	# Saves the chiptype into object lib
	lib<-phenodata$chiptype[1]
	lib<-as.character(lib)
}

# Account for the fact that annotation packages are from version 2.3 of Bioconductor
# named with an ".db" suffix. Add the suffix when missing to support data files
# from Chipster 1.3 and earlier. 
if (length(grep(".db", lib)) == 0 & length(grep("pmcdf", lib)) == 0) {
	lib <- paste(lib, ".db", sep="")
}

# Loads the correct annotation library
library(package=lib, character.only=T)

# Loads the data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Experimental design
groups<-phenodata[,pmatch(column,colnames(phenodata))]

lib2<-sub('.db','',lib)
env<-paste(lib2, "PATH", sep="") #KEGG identifiers
#env<-paste(lib2, "PFAM", sep="") #PFAM identifiers
#env<-paste(lib2, "GO2ALLPROBES", sep="") #GO identifiers
alleg<-get(env)
alleg<-as.list(alleg)
cmatrix<-getCmatrix(gene.list=alleg, present.genes=rownames(dat), min.size=minimum.category.size)
cmatrix$col.names <- paste("KEGG:", cmatrix$col.names, sep = "")

#cmatrix<-getCmatrix(keyword.list= alleg, present.genes = rownames(dat), min.size = 10, GO.ont = "CC") #GO.ont="BP", GO.ont="MF", GO.ont="CC"

#if(which.ontology=="GO") {
#	lib2<-sub('.db','',lib)
#	env<-paste(lib2, "GO2ALLPROBES", sep="")
#	alleg<-get(env)
#	alleg<-as.list(alleg)
#	tcmatrix<-getCmatrix(alleg)
#	tcmatrix<-tcmatrix[,dimnames(tcmatrix)[[2]] %in% rownames(dat2)]
#	cmatrix<-matrix(0, dim(dat2)[[1]], dim(tcmatrix)[[1]])
#	dimnames(cmatrix)<-list(dimnames(dat2)[[1]], dimnames(tcmatrix)[[1]]) 
#	cmatrix[match(dimnames(tcmatrix)[[2]], dimnames(dat2)[[1]]), ] <- t(tcmatrix)
#	cmatrix2<-cmatrix[,apply(cmatrix, 2, sum) >= minimum.category.size]
#}

# Runs the analysis
if(num.or.cat == "factor" & length(unique(groups)) == 2) {
	result<-safe(dat2, groups, cmatrix, alpha=p.value.threshold, min.size=minimum.category.size, local="t.Student")
} else if (num.or.cat == "factor" & length(unique(groups)) > 2) {
	result<-safe(dat2, groups, cmatrix, alpha=p.value.threshold, min.size=minimum.category.size, local="f.ANOVA")
} else if (num.or.cat == "continuous" & length(unique(groups)) >= 2) {
	for(i in 1:length(unique(groups))) {
		if(is.numeric(unique(groups)[i])==FALSE) {
			stop("CHIPSTER-NOTE: Your phenodata type continuous, but one of your group variables is not a number")
		}	
	}
	result<-safe(as.matrix(dat2), groups, cmatrix, alpha=p.value.threshold, min.size=minimum.category.size, local="t.LM")
} else {
	stop("CHIPSTER-NOTE: You need to have at least two groups to run this analysis")
}

if(length(which(result@global.pval <= p.value.threshold))==0) {
	bitmap(file="safeplot.png", width=w/72, height=h/72)
	plot(1, 1, col=0)
	text(1, 1, "This is a dummy image.", col=1)
	text(1, 0.9, "This has been generated, because no significant results were found.", col=1)
	text(1, 0.8, "These things happen.", col=1)
	dev.off()
	result.table <- data.frame(category=names(result@global.pval), size=rowSums(as.matrix(t(result@C.mat))), p.value=result@global.pval, p.adjusted=round(result@global.error, digits=3))
	result.table <- result.table[result.table$p.value<p.value.threshold,]
	result.table <- result.table[order(result.table$p.value),]
	write.table(result.table, file="safe.tsv", sep="\t", row.names=F, col.names=T, quote=F)
} else {
	# Plotting the results
	bitmap(file="safeplot.png", width=w/72, height=h/72)
	#Extreme=false means that all gene names are printed
	safeplot(result, extreme=FALSE)
	dev.off()
	# Writing a result table
	result.table <- data.frame(category=names(result@global.pval), size=rowSums(as.matrix(t(result@C.mat))), p.value=result@global.pval, p.adjusted=round(result@global.error, digits=3))
	result.table <- result.table[result.table$p.value<p.value.threshold,]
	result.table <- result.table[order(result.table$p.value),]
	write.table(result.table, file="safe.tsv", sep="\t", row.names=F, col.names=T, quote=F)
}
