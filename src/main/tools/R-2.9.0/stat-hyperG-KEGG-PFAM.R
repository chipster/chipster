# ANALYSIS Pathways/"Hypergeometric test for KEGG or PFAM" (Finds KEGG pathways or PFAM motifs that are over- or 
# under-represented.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT hypergeo.html
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value threshold)
# PARAMETER over.or.under.representation [over, under] DEFAULT over (Should over- or under-represented classes be seeked?)
# PARAMETER which.ontology [KEGG, PFAM] DEFAULT KEGG (Which ontology to use in the test?)

# Hypergeometrix test of gene enrichment to term categories
# Dario Greco 7.1.2007
# JTT 30.7.2007 (with heavy modifications) 

# Loads the libraries
library(GOstats)
library(KEGG.db)
library(PFAM.db)

# Renaming variables
pcut<-p.value.threshold
choisedirec<-over.or.under.representation
ontology<-which.ontology

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
if(ontology=="KEGG") {
   params<-new("KEGGHyperGParams", geneIds=myids, annotation=annotpkg, pvalueCutoff=pcut, testDirection=choisedirec)
   result<-hyperGTest(params)
}

if(ontology=="PFAM") {
   params<-new("PFAMHyperGParams", geneIds=myids, annotation=annotpkg, pvalueCutoff=pcut, testDirection=choisedirec)
   result<-hyperGTest(params)
}

# Writes the results to an HTML file
if(ontology=="KEGG" & sum(pvalues(result)<pcut)>=1) {
   ID<-names(pvalues(result)[pvalues(result)<=pcut])
   Pvalue<-pvalues(result)[pvalues(result)<=pcut]
   OddsRatio<-oddsRatios(result)[pvalues(result)<=pcut]
   ExpCount<-expectedCounts(result)[pvalues(result)<=pcut]
   Count<-geneCounts(result)[pvalues(result)<=pcut]
   Size<-universeCounts(result)[pvalues(result)<=pcut]
   Term<-as.vector(t(as.data.frame(mget(names(result@oddsRatios[result@pvalues<=pcut]), KEGGPATHID2NAME)))[,1])   
   write(x="<HTML>", file="hypergeo.html", append=T) 
   write(x="<BODY>", file="hypergeo.html", append=T)
   write(x="<TABLE border=1>", file="hypergeo.html", append=T)
   write(x="<CAPTION> KEGG - test for over-representation </CAPTION>", file="hypergeo.html", append=T)
   write(x="<TR> <TH>  </TH> <TH> ID </TH> <TH> Pvalue </TH> <TH> OddsRatio </TH> <TH> ExpCount </TH> <TH> Count </TH> <TH> Size </TH> <TH> Term </TH>  </TR>", file="hypergeo.html", append=T)
   for(i in 1:length(ID)) {
      write(x=paste("<TR> <TD> ", i, " </TD> <TD> ", ID[i], " </TD> <TD> ", round(Pvalue[i], digits=6), " </TD> <TD> ", round(OddsRatio[i], digits=2), " </TD> <TD> ", round(ExpCount[i], digits=2), " </TD> <TD> ", Count[i], " </TD> <TD> ", Size[i], " </TD> <TD> ", "<a href=", cat("\""), "http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=gn&dbkey=pathway&keywords=", ID[i], "&page1", cat("\""), ">", Term[i], "</a>", "</TD> </TR>", sep=""), file="hypergeo.html", append=T)
   }
   write(x="</BODY>", file="hypergeo.html", append=T)
   write(x="</HTML>", file="hypergeo.html", append=T)      
} 
if(ontology=="KEGG" & sum(pvalues(result)<pcut)<1) {
   write(x="<HTML>", file="hypergeo.html", append=T) 
   write(x="<BODY>", file="hypergeo.html", append=T)
   write(x="No significant results found! <br>", file="hypergeo.html", append=T)
   write(x="</BODY>", file="hypergeo.html", append=T)
   write(x="</HTML>", file="hypergeo.html", append=T)      
}

if(ontology=="PFAM" & sum(pvalues(result)<pcut)>=1) {
   ID<-names(pvalues(result)[pvalues(result)<=pcut])
   Pvalue<-pvalues(result)[pvalues(result)<=pcut]
   OddsRatio<-oddsRatios(result)[pvalues(result)<=pcut]
   ExpCount<-expectedCounts(result)[pvalues(result)<=pcut]
   Count<-geneCounts(result)[pvalues(result)<=pcut]
   Size<-universeCounts(result)[pvalues(result)<=pcut]
   write(x="<HTML>", file="hypergeo.html", append=T) 
   write(x="<BODY>", file="hypergeo.html", append=T)
   write(x="<TABLE border=1>", file="hypergeo.html", append=T)
   write(x="<CAPTION> PFAM - test for over-representation </CAPTION>", file="hypergeo.html", append=T)
   write(x="<TR> <TH>  </TH> <TH> ID </TH> <TH> Pvalue </TH> <TH> OddsRatio </TH> <TH> ExpCount </TH> <TH> Count </TH> <TH> Size </TH> <TH> Term </TH>  </TR>", file="hypergeo.html", append=T)
   for(i in 1:length(ID)) {
      write(x=paste("<TR> <TD> ", i, " </TD> <TD> ", ID[i], " </TD> <TD> ", round(Pvalue[i], digits=2), " </TD> <TD> ", round(OddsRatio[i], digits=2), " </TD> <TD> ", round(ExpCount[i], digits=2), " </TD> <TD> ", Count[i], " </TD> <TD> ", Size[i], " </TD> <TD> ", "<a href=", cat("\""), "http://www.sanger.ac.uk/cgi-bin/Pfam/qquerypfam.pl?terms=", ID[i], cat("\""), ">", ID[i], "</a>", "</TD> </TR>", sep=""), file="hypergeo.html", append=T)
   }
   write(x="</BODY>", file="hypergeo.html", append=T)
   write(x="</HTML>", file="hypergeo.html", append=T)      
} 
if(ontology=="PFAM" & sum(pvalues(result)<pcut)<1) {
   write(x="<HTML>", file="hypergeo.html", append=T) 
   write(x="<BODY>", file="hypergeo.html", append=T)
   write(x="No significant results found! <br>", file="hypergeo.html", append=T)
   write(x="</BODY>", file="hypergeo.html", append=T)
   write(x="</HTML>", file="hypergeo.html", append=T)      
}

