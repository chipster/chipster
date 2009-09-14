# ANALYSIS Annotation/"Agilent, Affymetrix or Illumina genelist" (Annotates the genes, and creates a web page of the result. 
# Currently, this function only works with Agilent, Affymetrix and Illumina data if annotation packages for the
# used chip types are available.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT annotations.html, annotations.tsv
# PARAMETER expression.column COLUMN_SEL DEFAULT EMPTY (Column that contains the expression values)
# PARAMETER p.value.column COLUMN_SEL DEFAULT EMPTY (Column that contains the p-values)


# Creates an HTML file of Affymetrix genelist annotations
# JTT 13.6.2006

# Reads the chiptype from phenodata table
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
if(phenodata$chiptype[1]!="cDNA" | phenodata$chiptype[1]!="Illumina") {
   # Saves the chiptype into object lib
   lib<-phenodata$chiptype[1]
   lib<-as.character(lib)
}
lib <- paste(lib, ".db", sep="")

# Loads the correct annotation library
library(package=lib, character.only=T)
library(annaffy)

# Loads the data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Creating annotations from the library
if(expression.column=="EMPTY" & p.value.column=="EMPTY") {
   annot.cols<-aaf.handler()
   annot.table<-aafTableAnn(row.names(dat), lib, annot.cols)
   saveHTML(annot.table, "annotations.html", title="Annotations for the selected gene list")
   saveText(annot.table, "annotations.tsv")
}

if(expression.column!="EMPTY" & p.value.column=="EMPTY") {
   annot.cols<-aaf.handler()
   annot.table<-aafTableAnn(row.names(dat), lib, annot.cols)
   expression<-data.frame(FoldChange=dat[,grep(expression.column, colnames(dat))])
   rownames(expression)<-rownames(dat)
   expression<-aafTable(items=expression, signed=T)
   expression@probeids<-rownames(dat)
   annot.table2<-merge(annot.table, expression)
   saveHTML(annot.table2, "annotations.html", title="Annotations for the selected gene list")
   saveText(annot.table2, "annotations.tsv")
}

if(expression.column=="EMPTY" & p.value.column!="EMPTY") {
   annot.cols<-aaf.handler()
   annot.table<-aafTableAnn(row.names(dat), lib, annot.cols)
   pvalues<-data.frame(P.Adjusted=dat[,grep(p.value.column, colnames(dat))])
   rownames(pvalues)<-rownames(dat)
   pvalues<-aafTable(items=pvalues)
   pvalues@probeids<-rownames(dat)
   annot.table2<-merge(annot.table, pvalues)
   saveHTML(annot.table2, "annotations.html", title="Annotations for the selected gene list")
   saveText(annot.table2, "annotations.tsv")
}

if(expression.column!="EMPTY" & p.value.column!="EMPTY") {
   annot.cols<-aaf.handler()
   annot.table<-aafTableAnn(row.names(dat), lib, annot.cols)
   expression<-data.frame(FoldChange=dat[,grep(expression.column, colnames(dat))])
   rownames(expression)<-rownames(dat)
   expression<-aafTable(items=expression, signed=T)
   expression@probeids<-rownames(dat)
   annot.table2<-merge(annot.table, expression)
   pvalues<-data.frame(P.Adjusted=dat[,grep(p.value.column, colnames(dat))])
   rownames(pvalues)<-rownames(dat)
   pvalues<-aafTable(items=pvalues)
   pvalues@probeids<-rownames(dat)
   annot.table3<-merge(annot.table2, pvalues)
   saveHTML(annot.table3, "annotations.html", title="Annotations for the selected gene list")
   saveText(annot.table3, "annotations.tsv")
}


