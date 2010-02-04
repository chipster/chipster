# ANALYSIS Pathways/"Hypergeometric test for GO" (Finds GO ontology classes that are over- or under-represented.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT hypergeo.html
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value threshold)
# PARAMETER over.or.under.representation [over, under] DEFAULT over (Should over- or under-represented classes be seeked?)

# Hypergeometrix test of gene enrichment to term categories
# Dario Greco 7.1.2007
# JTT 30.7.2007 (with heavy modifications) 

# Loads the libraries
library(GOstats)

# Renaming variables
pcut<-p.value.threshold
choisedirec<-over.or.under.representation

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Reading phenodata for chiptype
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

# Loads the correct Affymetrix annotation library
library(package=lib, character.only=T)

# in the next 4 rows it converts the microarray IDs into Entrez Gene IDs
lib2<-sub('.db','',lib)
env<-paste(lib2, "ENTREZID", sep="")
alleg<-get(env)
alleg<-as.list(alleg)
alleg<-unlist(alleg)
alleg<-as.data.frame(alleg)
myids<-unique(alleg[rownames(dat),])

#CHOOSE THE ANNOTATION PACKAGE
annotpkg<-lib

#Account for the fact that updated alternative CDF annotation packages for Affymetrix
#hgu133A, hgu133A2 are no longer supported and should be replaced with hgu133plus2 
if (annotpkg=="hgu133ahsentrezg.db" || "hgu133a2hsentrezg.db") {
	annotpkg <- "hgu133plus2hsentrezg.db"
}

params<-new("GOHyperGParams", geneIds=myids, annotation=annotpkg, ontology="BP", pvalueCutoff=pcut, conditional=TRUE,testDirection=choisedirec)
resultBP<-hyperGTest(params)
params<-new("GOHyperGParams", geneIds=myids, annotation=annotpkg, ontology="MF", pvalueCutoff=pcut, conditional=TRUE,testDirection=choisedirec)
resultMF<-hyperGTest(params)
params<-new("GOHyperGParams", geneIds=myids, annotation=annotpkg, ontology="CC", pvalueCutoff=pcut, conditional=TRUE,testDirection=choisedirec)
resultCC<-hyperGTest(params)

if( sum(pvalues(resultCC)<pcut) + sum(pvalues(resultBP)<pcut) + sum(pvalues(resultMF)<pcut) >=1) {
   htmlReport(resultBP, "hypergeo.html", append=T)
   htmlReport(resultMF, "hypergeo.html", append=T)
   htmlReport(resultCC, "hypergeo.html", append=T)
}

if( sum(pvalues(resultCC)<pcut) + sum(pvalues(resultBP)<pcut) + sum(pvalues(resultMF)<pcut) <1) {
   write(x="<HTML>", file="hypergeo.html", append=T) 
   write(x="<BODY>", file="hypergeo.html", append=T)
   write(x="No significant results found! <br>", file="hypergeo.html", append=T)
   write(x="</BODY>", file="hypergeo.html", append=T)
   write(x="</HTML>", file="hypergeo.html", append=T)   
}
