# ANALYSIS Pathways/"Hypergeometric test for GO" (Finds GO ontology classes that are over- or under-represented.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv
# OUTPUT hypergeo.tsv, hypergeo.html
# PARAMETER ontology [all, biological_process, molecular_function, cellular_component] DEFAULT biological_process (the ontology to be analyzed)
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value threshold)
# PARAMETER over.or.under.representation [over, under] DEFAULT over (Should over or under-represented classes be seeked?)

# Hypergeometrix test of gene enrichment to term categories
# Dario Greco 7.1.2007
# JTT 30.7.2007 (with heavy modifications) 
#
# modified MG 5.2.2010
# modified IS 14.9.2010

# Loads the libraries
library(GOstats)

# Loads the normalized data
dat<-read.table("normalized.tsv", header=TRUE, sep="\t", row.names=1)

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

#Account for the fact that updated alternative CDF annotation packages for some Affymetrix
#arrays are no longer supported and should be replaced with their Bioconductor counterpart 
if (annotpkg=="hgu133a2hsentrezg.db") {
	annotpkg <- "hgu133a2.db"
}
if (annotpkg=="hgu133ahsentrezg.db") {
	annotpkg <- "hgu133a.db"
}
if (annotpkg=="ath1121501hsentrezg.db") {
	annotpkg <- "ath1121501.db"
}
if (annotpkg=="aghsentrezg.db") {
	annotpkg <- "ag.db"
}
if (annotpkg=="ygs98scentrezg.db") {
	annotpkg <- "ygs98.db"
}
if (annotpkg=="yeast2scentrezg.db") {
	annotpkg <- "yeast2.db"
}

# define the output variable
output <- data.frame(total=integer(0), expectation=numeric(0), observation=integer(0), p.value=numeric(0), description=character(0), ontology=character(0))

if (ontology == 'biological_process' || ontology == 'all') {
  params<-new("GOHyperGParams", geneIds=myids, annotation=annotpkg, ontology="BP", pvalueCutoff=p.value.threshold, conditional=TRUE,testDirection=over.or.under.representation)
  go <- hyperGTest(params)
  go.table <- summary(go)
  if (nrow(go.table)>0) {
    rownames(go.table) <- go.table[,1]
    go.table <- go.table[,c(6, 4, 5, 2, 7)]
    go.table$ontology <- 'biological process'
    colnames(go.table) <- colnames(output)
    output <- rbind(output, go.table)
    htmlReport(go, file='hypergeo.html', append=TRUE)
  }
}

if (ontology == 'molecular_function' || ontology == 'all') {
  params<-new("GOHyperGParams", geneIds=myids, annotation=annotpkg, ontology="MF", pvalueCutoff=p.value.threshold, conditional=TRUE,testDirection=over.or.under.representation)
  go <- hyperGTest(params)
  go.table <- summary(go)
  if (nrow(go.table)>0) {
    rownames(go.table) <- go.table[,1]
    go.table <- go.table[,c(6, 4, 5, 2, 7)]
    go.table$ontology <- 'molecular function'
    colnames(go.table) <- colnames(output)
    output <- rbind(output, go.table)
    htmlReport(go, file='hypergeo.html', append=TRUE)
  }
}

if (ontology == 'cellular_component' || ontology == 'all') {
  params<-new("GOHyperGParams", geneIds=myids, annotation=annotpkg, ontology="CC", pvalueCutoff=p.value.threshold, conditional=TRUE,testDirection=over.or.under.representation)
  go <- hyperGTest(params)
  go.table <- summary(go)
  if (nrow(go.table)>0) {
    rownames(go.table) <- go.table[,1]
    go.table <- go.table[,c(6, 4, 5, 2, 7)]
    go.table$ontology <- 'cellular component'
    colnames(go.table) <- colnames(output)
    output <- rbind(output, go.table)
    htmlReport(go, file='hypergeo.html', append=TRUE)
  }
}

# write outputs
write.table(output, file='hypergeo.tsv', quote=FALSE, sep='\t')
if (nrow(output)==0)
  write('<html>\n\t<body>\n\t\tNo significant results found!</br />\n\t</body>\n</html>', file='hypergeo.html')

# EOF