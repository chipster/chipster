# ANALYSIS Pathways/"Gene set test" (Testing the significance of a bunch of genes as group, not separately.
# This tool makes a global test of gene expression, and gives the significance of the expression 
# of a group of genes. You can use your own genelist or search for the most sigficant KEGG or GO pathways. Searching
# for the most significant pathways works only, if the whole dataset is used, i.e, data should not be prefiltered.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT multtest.png, globaltest-result-table.tsv
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER pathway.or.genelist [KEGG, GO, current] DEFAULT current (Which gene list to test)
# PARAMETER use.multiple.testing.correction [yes, no] DEFAULT yes (Should multiple testing correction be used)
# PARAMETER number.of.groups.to.visualize [4, 9, 16] DEFAULT 4 (Number of most significant groups to visualize)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Testing the statistical significance of several genes as a group or using KEGG or GO pathways
# This only works with KEGG and GO for the whole normalized data (no pre-filtering allowed)
# JTT 14.6.2006

# Loading the libraries
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
if(phenodata$chiptype[1]!="cDNA" | phenodata$chiptype[1]!="Illumina") {
   # Saves the chiptype into object lib
   lib<-phenodata$chiptype[1]
   lib<-as.character(lib)
}

library(package=lib, character.only=T)
library(KEGG)
library(GO)
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
groups<-phenodata[,grep(column, colnames(phenodata))]

# Actual testing and plotting
x<-as.numeric(x) 

# Testing with KEGG pathways
if(pathways=="KEGG") {
   # Calculating the test
   env<-paste(lib, "PATH2PROBE", sep="")
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
   table.out<-data.frame(pathway=names(test.kegg), pvalue=p.value(test.kegg))
   names(test.kegg)<-as.list(KEGGPATHID2NAME)[names(test.kegg)]
   table.out<-data.frame(table.out, Description=names(test.kegg))
   table.out<-table.out[order(table.out$pvalue),]
   write.table(table.out, file="globaltest-result-table.tsv", sep="\t", row.names=T, col.names=T, quote=F)
   s<-c(sqrt(x), sqrt(x))
   bitmap(file="multtest.png", width=w/72, height=h/72)
   split.screen(c(sqrt(x), sqrt(x)))
   for(i in 1:x) {
      y<-geneplot(test.kegg[i], plot=F, addlegend=F)
      env2<-paste(lib, "SYMBOL", sep="")
      env2<-get(env2)
      names(y)<-as.list(env2)[names(y)]
      screen(i)
      par(mar=c(2,0,2,0)+0.1, yaxt="n")
      plot(y, main=substr(names(test.kegg[i]), 1, 25), cex.main=0.75)
  }
  dev.off()
}

# Testing with GO pathways
if(pathways=="GO") {
   # The same with GO-terms
   env<-paste(lib, "GO2ALLPROBES", sep="")
   go2allprobes<-get(env)
   go<-as.list(go2allprobes)
   test.go<-globaltest(as.matrix(dat2), groups, go)
   test.go<-sort(test.go)
   test.go2<-test.go[1:x,]
   table.out<-data.frame(pathway=names(test.go2), pvalue=p.value(test.go2))
   n<-c()
   for(i in 1:x) {
      n<-c(n, get(names(test.go2[i,]),GOTERM)@Term)
   }
   names(test.go2)<-n
   table.out<-data.frame(table.out, Description=names(test.go2))
   table.out<-table.out[order(table.out$pvalue),]
   write.table(table.out, file="globaltest-result-table.tsv", sep="\t", row.names=T, col.names=T, quote=F)
   s<-c(sqrt(x), sqrt(x)) 
   bitmap(file="multtest.png", width=w/72, height=h/72)
   split.screen(c(sqrt(x), sqrt(x)))
   for(i in 1:x) {
      y<-geneplot(test.go[1], plot=F, addlegend=F)
      env2<-paste(lib, "SYMBOL", sep="")
      env2<-get(env2)
      names(y)<-as.list(env2)[names(y)]
      screen(i)
      par(mar=c(2,0,2,0)+0.1, yaxt="n")
      plot(y, main=substr(names(test.go2[i]), 1, 25), cex.main=0.75)
   }
   dev.off()
}

if(pathways=="current") {
   test.current<-globaltest(X=as.matrix(dat2), Y=groups, genesets=rep(1, nrow(dat2)))
   bitmap(file="multtest.png", width=w/72, height=h/72)
   geneplot(test.current)
   dev.off()
   write(x="Dont't worry! You've run the analysis using current gene list. Thus, this file should be empty.", file="globaltest-result-table.tsv")
}
