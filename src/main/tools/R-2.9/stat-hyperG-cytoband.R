# ANALYSIS Pathways/"Hypergeometric test for cytobands" (Finds chromosomal locations, i.e., cytobands for which 
# the genes in the selected dataset are either over- or under-represented.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT hypergeo.html
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value threshold)
# PARAMETER over.or.under.representation [over, under] DEFAULT over (Should over- or under-represented classes be seeked?)

# Hypergeometrix test of gene enrichment to term categories
# JTT 22.1.2009

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
if (length(grep(".db", lib)) == 0) {
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

# Runs the tests
params<-new("ChrMapHyperGParams", geneIds=myids, annotation=annotpkg, pvalueCutoff=pcut, testDirection=choisedirec, conditional=FALSE)   
result<-hyperGTest(params)

# Writes the results to an HTML file
if(sum(pvalues(result)<pcut)>=1) {
   ID<-names(pvalues(result)[pvalues(result)<=pcut])
   Pvalue<-pvalues(result)[pvalues(result)<=pcut]
   OddsRatio<-oddsRatios(result)[pvalues(result)<=pcut]
   ExpCount<-expectedCounts(result)[pvalues(result)<=pcut]
   Count<-geneCounts(result)[pvalues(result)<=pcut]
   Size<-universeCounts(result)[pvalues(result)<=pcut]
   write(x="<HTML>", file="hypergeo.html", append=T) 
   write(x="<BODY>", file="hypergeo.html", append=T)
   write(x="<TABLE border=1>", file="hypergeo.html", append=T)
   write(x="<CAPTION> Cytoband - test for over-representation </CAPTION>", file="hypergeo.html", append=T)
   write(x="<TR> <TH>  </TH> <TH> Cytoband </TH> <TH> Pvalue </TH> <TH> OddsRatio </TH> <TH> ExpCount </TH> <TH> Count </TH> <TH> Size </TH>  </TR>", file="hypergeo.html", append=T)
   for(i in 1:length(ID)) {
      write(x=paste("<TR> <TD> ", i, " </TD> <TD> ", ID[i], " </TD> <TD> ", round(Pvalue[i], digits=6), " </TD> <TD> ", round(OddsRatio[i], digits=2), " </TD> <TD> ", round(ExpCount[i], digits=2), " </TD> <TD> ", Count[i], " </TD> <TD> ", Size[i], " </TD> </TR>", sep=""), file="hypergeo.html", append=T)
   }
   write(x="</BODY>", file="hypergeo.html", append=T)
   write(x="</HTML>", file="hypergeo.html", append=T)      
} 
if(sum(pvalues(result)<pcut)<1) {
   write(x="<HTML>", file="hypergeo.html", append=T) 
   write(x="<BODY>", file="hypergeo.html", append=T)
   write(x="No significant results found! <br>", file="hypergeo.html", append=T)
   write(x="</BODY>", file="hypergeo.html", append=T)
   write(x="</HTML>", file="hypergeo.html", append=T)      
}

