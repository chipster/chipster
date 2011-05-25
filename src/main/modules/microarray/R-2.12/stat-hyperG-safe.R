# TOOL stat-hyperG-safe.R: "SAFE test for KEGG pathway enrichment" (Finds KEGG pathways that are over- or under-represented in the selected gene list.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT safeplot.png: safeplot.png 
# OUTPUT safe.tsv: safe.tsv 
# PARAMETER column: column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER p.value.threshold: p.value.threshold TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)
# PARAMETER which.ontology: which.ontology TYPE [KEGG: KEGG, GO: GO] DEFAULT KEGG (Which ontology to use in the test?)
# PARAMETER minimum.category.size: minimum.category.size TYPE INTEGER FROM 1 TO 100 DEFAULT 10 (Minimum size for categories to be evaluated)

# Modified 5.11.2009, MG
# Accounting for changes in SparseM package regarding matrix.csr class

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

# Creates a C matrix
if(which.ontology=="KEGG") {
	lib2<-sub('.db','',lib)
	env<-paste(lib2, "PATH", sep="")
	alleg<-get(env)
	alleg<-as.list(alleg)
	cmatrix<-getCmatrix(gene.list=alleg, present.genes=rownames(dat), min.size=minimum.category.size)
	#cmatrix<-getCmatrix(gene.list=alleg, present.genes=rownames(dat))
}

if(which.ontology=="GO") {
	env<-paste(lib, "GO2ALLPROBES", sep="")
	alleg<-get(env)
	alleg<-as.list(alleg)
	tcmatrix<-getCmatrix(alleg)
	tcmatrix<-tcmatrix[,dimnames(tcmatrix)[[2]] %in% rownames(dat2)]
	cmatrix<-matrix(0, dim(dat2)[[1]], dim(tcmatrix)[[1]])
	dimnames(cmatrix)<-list(dimnames(dat2)[[1]], dimnames(tcmatrix)[[1]]) 
	cmatrix[match(dimnames(tcmatrix)[[2]], dimnames(dat2)[[1]]), ] <- t(tcmatrix)
	cmatrix2<-cmatrix[,apply(cmatrix, 2, sum) >= minimum.category.size]
}

# Runs the analysis
result<-safe(dat2, groups, cmatrix, alpha=p.value.threshold, min.size=minimum.category.size)

# Are there any significant results?
# The following check is not valid any more: you always get p-values, but they are just below the threshold. 
# Instead one should check if result@C.mat exists?
if(length(result@global.pval)==0) {
	bitmap(file="safeplot.png", width=w/72, height=h/72)
	plot(1, 1, col=0)
	text(1, 1, "This is a dummy image.", col=1)
	text(1, 0.9, "This has been generated, because no significant results were found.", col=1)
	text(1, 0.8, "These things happen.", col=1)
	dev.off()
	result.table <- data.frame(Category=names(result@global.pval), Size=rowSums(as.matrix(t(result@C.mat))), P.Value=result@global.pval)
	result.table <- result.table[result.table$P.Value<p.value.threshold,]
	result.table <- result.table[order(result.table$P.Value),]
	write.table(result.table, file="safe.tsv", sep="\t", row.names=F, col.names=T, quote=F)
} else {
	# Plotting the results
	bitmap(file="safeplot.png", width=w/72, height=h/72)
	safeplot(result)
	dev.off()
	# Writing a result table
	result.table <- data.frame(Category=names(result@global.pval), Size=rowSums(as.matrix(t(result@C.mat))), P.Value=result@global.pval)
	result.table <- result.table[result.table$P.Value<p.value.threshold,]
	result.table <- result.table[order(result.table$P.Value),]
	write.table(result.table, file="safe.tsv", sep="\t", row.names=F, col.names=T, quote=F)
}

