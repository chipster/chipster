# TOOL stat-geneset.R: "Gene set test" (Testing the significance of a bunch of genes as group, not separately. This tool makes a global test of gene expression, and gives the significance of the expression of a group of genes. You can use your own genelist or search for the most sigficant KEGG or GO pathways. Searching for the most significant pathways works only, if the whole dataset is used, i.e, data should not be prefiltered.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT multtest.pdf: multtest.pdf 
# OUTPUT globaltest-result-table.tsv: globaltest-result-table.tsv 
# PARAMETER column: column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER pathway.or.genelist: pathway.or.genelist TYPE [KEGG: KEGG, GO: GO, current: current] DEFAULT KEGG (Which gene list to test)
# PARAMETER use.multiple.testing.correction: use.multiple.testing.correction TYPE [yes: yes, no: no] DEFAULT yes (Should multiple testing correction be used)
# PARAMETER number.of.groups.to.visualize: number.of.groups.to.visualize TYPE [4: 4, 9: 9, 16: 16] DEFAULT 4 (Number of most significant groups to visualize)
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Testing the statistical significance of several genes as a group or using KEGG or GO pathways
# This only works with KEGG and GO for the whole normalized data (no pre-filtering allowed)
# JTT 14.6.2006
#
# modified by MG 31.12.2009
#
# modified by MG, 3.5.2010 to exclude single gene genesets and added group labels in plots
# added more info columns in results table

#column<-"group"
#pathway.or.genelist <-"KEGG"
#use.multiple.testing.correction <-"yes"
#number.of.groups.to.visualize <-16
#image.width <-600
#image.height <-600

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
mt<-use.multiple.testing.correction
x<-number.of.groups.to.visualize
w<-image.width
h<-image.height

# Reading in data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Needs a vector groups that specifies which sample to compare
groups<-phenodata[,pmatch(column,colnames(phenodata))]

# Actual testing and plotting
x<-as.numeric(x) 

# Testing with KEGG pathways
if(pathways=="KEGG") {
   # Calculating the test
   lib2<-sub('.db','',lib)
   env<-paste(lib2, "PATH2PROBE", sep="")
   pathway2probe<-get(env)
   kegg<-as.list(pathway2probe)
   test.kegg<-globaltest(as.matrix(dat2), groups, kegg)
   if(mt=="yes") {
      test.kegg<-gt.multtest(test.kegg)
      test.kegg<-sort(test.kegg)
   }
   if(mt=="no") {
      test.kegg<-sort(test.kegg)
   }
#   table.out<-data.frame(pathway=names(test.kegg), pvalue=p.value(test.kegg))
   table.out<-data.frame(pathwayID=names(test.kegg), test.kegg@res)
   names(test.kegg)<-as.list(KEGGPATHID2NAME)[names(test.kegg)]
   table.out<-data.frame(table.out, Description=names(test.kegg))
   table.out<-table.out[order(table.out$P.value),]
   # Find out if there are categories including a single gene
	# and if so remove them
   indices <- 1:length(table.out$Tested)
   indices <- indices[table.out$Tested>1]
   table.out <- table.out[table.out$Tested>1,]
   write.table(table.out, file="globaltest-result-table.tsv", sep="\t", row.names=T, col.names=T, quote=F)
   s<-c(sqrt(x), sqrt(x))
   pdf(file="multtest.pdf", width=w/72, height=h/72)
   split.screen(c(sqrt(x), sqrt(x)))
   if (length (indices) < x) {
	   x <- length(indices)
   } 
   for(i in 1:x) {
	   y<-geneplot(test.kegg[indices[i]], plot=F, addlegend=T)
	   env2<-paste(lib2, "SYMBOL", sep="")
	   env2<-get(env2)
	   names(y)<-as.list(env2)[names(y)]
	   screen(i)
	   par(mar=c(2,0,2,0)+0.1, yaxt="n", cex=0.75)
	   plot(y, main=substr(names(test.kegg[i]), 1, 25), cex.main=1.5, cex.axis=1.25)
   } 
  dev.off()
}

# Testing with GO pathways
if(pathways=="GO") {
	# The same with GO-terms
	lib2<-sub('.db','',lib)
	env<-paste(lib2, "GO2ALLPROBES", sep="")
	go2allprobes<-get(env)
	go<-as.list(go2allprobes)
	test.go<-globaltest(as.matrix(dat2), groups, go)
	
	if (mt=="yes") {
		test.go <- gt.multtest(test.go)
		test.go <- sort(test.go)
	}
	if (mt=="no") {
		test.go <- sort(test.go)
	}

	
#	test.go<-sort(test.go)
	test.go2<-test.go[1:x,]
	table.out<-data.frame(ontologyID=names(test.go2), test.go2@res)
	n<-c()
	for(i in 1:x) {
		n<-c(n, get(names(test.go2[i,]),GOTERM)@Term)
	}
	names(test.go2)<-n
	table.out<-data.frame(table.out, Description=names(test.go2))
	table.out<-table.out[order(table.out$P.value),]
# Find out if there are categories including a single gene
# and if so remove them
	indices <- 1:length(table.out$Tested)
	indices <- indices[table.out$Tested>1]
	table.out <- table.out[table.out$Tested>1,]
	write.table(table.out, file="globaltest-result-table.tsv", sep="\t", row.names=T, col.names=T, quote=F)
	s<-c(sqrt(x), sqrt(x)) 
	pdf(file="multtest.pdf", width=w/72, height=h/72)
	split.screen(c(sqrt(x), sqrt(x)))
	if (length (indices) < x) {
		x <- length(indices)
	} 
	for(i in 1:x) {
		y<-geneplot(test.go[indices[i]], plot=F, addlegend=T)
		env2<-paste(lib2, "SYMBOL", sep="")
		env2<-get(env2)
		names(y)<-as.list(env2)[names(y)]
		screen(i)
		par(mar=c(2,0,2,0)+0.1, yaxt="n", cex=0.5)
		plot(y, main=substr(names(test.go2[i]), 1, 25), cex.main=1.5, cex.axis=1.25)
	}
	dev.off()
}

if(pathways=="current") {
   test.current<-globaltest(X=as.matrix(dat2), Y=groups, genesets=rep(1, nrow(dat2)))
   pdf(file="multtest.pdf", width=w/72, height=h/72)
   geneplot(test.current)
   dev.off()
   write(x="Dont't worry! You've run the analysis using current gene list. Thus, this file should be empty.", file="globaltest-result-table.tsv")
}
