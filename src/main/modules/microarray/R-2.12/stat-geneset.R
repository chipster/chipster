# TOOL stat-geneset.R: "Gene set test" (Testing the significance of a bunch of genes as group, not separately. This tool makes a global test of gene expression, and gives the significance of the expression of a group of genes. You can use your own genelist or search for the most sigficant KEGG or GO pathways. Searching for the most significant pathways works only, if the whole dataset is used, i.e, data should not be prefiltered.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT multtest.pdf: multtest.pdf 
# OUTPUT globaltest-result-table.tsv: globaltest-result-table.tsv 
# PARAMETER column: column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER pathway.or.genelist: "Pathway or genelist" TYPE [KEGG: KEGG, GO: GO, current: current] DEFAULT KEGG (Which gene list to test)
# PARAMETER mult.test.cor: "Multiple testing method" TYPE [Holm: Holm, BH: BH, BY:BY] DEFAULT BH (Type of multiple testing correction to be used)
# PARAMETER minimum.category.size: "Minimum category size" TYPE INTEGER FROM 1 TO 100 DEFAULT 5 (Minimum size for categories to be evaluated)
# PARAMETER p.value.threshold: "P-value threshold" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER number.of.groups.to.visualize: "Number of groups to visualize" TYPE INTEGER FROM 1 TO 100 DEFAULT 5 (Number of most significant groups to visualize)

# JTT 14.6.2006: Testing the statistical significance of several genes as a group or using KEGG or GO pathways. This only works with KEGG and GO for the whole normalized data (no pre-filtering allowed)
# MG 31.12.2009:
# MG 03.05.2010: to exclude single gene genesets and added group labels in plots. Added more info columns in results table
# MK 03.09.2013: updated to R.3.0 which does not support globaltest / geneplot functions

#column<-"group"
#pathway.or.genelist <-"KEGG"
#use.multiple.testing.correction <-"yes"
#number.of.groups.to.visualize <-16

# Loading the libraries
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

library(package=lib, character.only=T)
library(KEGG.db)
library(GO.db)
library(globaltest)

# Renaming variables
pathways<-pathway.or.genelist
x<-as.numeric(number.of.groups.to.visualize)

# Reading in data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[, grep("flag", names(dat))]
dat2<-dat[, grep("chip", names(dat))]

# Needs a vector groups that specifies which sample to compare
groups<-phenodata[, pmatch(column,colnames(phenodata))]


#mapping probe information to ENTREZIDs
#ls("package:hgug4851a.db")
lib2<-sub('.db','',lib)
env<-paste(lib2, "ENTREZID", sep="")
probe2entrez <- get(env)
env<-as.list(probe2entrez)

# Testing with KEGG pathways
if(pathways=="KEGG") {
	env2<-paste(lib2, "PATH2PROBE", sep="")
	pathway2probe <- get(env2)
	env2<-as.list(pathway2probe)
	path_count <- sapply(env2, length)
   
	gt.options(transpose=TRUE)
	test.gt <- gtKEGG(groups, as.matrix(dat2), probe2entrez = env, annotation="org.Hs.eg.db", multtest=mult.test.cor)

	if(mult.test.cor == "Holm") { mult.test.cor = "holm"; }
	table.out<-data.frame(pathwayID=names(test.gt), Genes=path_count[names(test.gt)], Tested=size(test.gt), 
			             Statistic.Q=test.gt@result[, "Statistic"], Expected.Q=test.gt@result[, "Expected"], sd.of.Q=test.gt@result[,"Std.dev"],
						 p.value=p.value(test.gt), p.adjusted=result(test.gt)[, mult.test.cor], Description=alias(test.gt))
	rownames(table.out)<-paste("KEGG:", rownames(table.out), sep = "")
}

# Testing with GO pathways
if(pathways=="GO") {
	env2<-paste(lib2, "GO2ALLPROBES", sep="")
	pathway2probe <- get(env2)
	env2<-as.list(pathway2probe)
	path_count <- sapply(env2, length)

	gt.options(transpose=TRUE)
	test.gt <- gtGO(groups, as.matrix(dat2), probe2entrez = env, annotation="org.Hs.eg.db", multtest=mult.test.cor)

	if(mult.test.cor == "Holm") { mult.test.cor = "holm"; }
	table.out<-data.frame(pathwayID=names(test.gt), Genes=path_count[names(test.gt)], Tested=size(test.gt), 
						  Statistic.Q=test.gt@result[, "Statistic"], Expected.Q=test.gt@result[, "Expected"], sd.of.Q=test.gt@result[,"Std.dev"],
						  p.value=p.value(test.gt), p.adjusted=result(test.gt)[, mult.test.cor], Description=alias(test.gt))
}

# Testing with the phenotype. Shows those genes that are differentially expressed between conditions
if(pathways=="current") {
	gt.options(transpose=TRUE)
	test.gt <- gt(groups, as.matrix(dat2))

	if(mult.test.cor == "Holm") { mult.test.cor = "holm"; }
	table.out<-data.frame(pathwayID=column, Genes=nrow(dat2), Tested=size(test.gt), 
		   Statistic.Q=test.gt@result[, "Statistic"], Expected.Q=test.gt@result[, "Expected"], sd.of.Q=test.gt@result[,"Std.dev"],
		   p.value=p.value(test.gt), p.adjusted=p.value(test.gt), Description="Phenodata column")
	rownames(table.out) <- column
}

# Make output and write
table.out <- table.out[order(table.out$p.value),]
table.out <- table.out[(which(table.out$Tested >= minimum.category.size)), ]
table.out <- table.out[(which(table.out$p.adjusted <= p.value.threshold)), ]
write.table(table.out, file="globaltest-result-table.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# Create plots
if(nrow(table.out) > 0) {
	if(pathways=="current") {
		ft <- features(test.gt, alias=paste(lib2, "SYMBOL", sep=""), pdf="multtest.pdf")
	} else {
		path.names <- as.vector(table.out[1:x, "pathwayID"])
  		ft <- features(test.gt[ path.names[!is.na(path.names)] ], alias=paste(lib2, "SYMBOL", sep=""), pdf="multtest.pdf")
	}
} else {
	pdf(file="multtest.pdf", width=600/72, height=600/72)
	plot(1, 1, col=0)
	text(1, 1, "This is a dummy image that was generated because no significant results were found", col=1)
	text(1, 0.8, "Please reset your p-value threshold to a larger value.", col=1)
	dev.off()	
}



