# TOOL stat-hyperG-GO.R: "Hypergeometric test for GO" (Finds GO ontology classes that are over- or under-represented.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT hypergeo.tsv: hypergeo.tsv 
# OUTPUT hypergeo.html: hypergeo.html 
# PARAMETER ontology: ontology TYPE [all: all, biological_process: biological_process, molecular_function: molecular_function, cellular_component: cellular_component] DEFAULT biological_process (The ontology to be analyzed.)
# PARAMETER p.value.threshold: p.value.threshold TYPE DECIMAL DEFAULT 0.05 (P-value threshold.)
# PARAMETER minimum.population: minimum.population TYPE INTEGER FROM 1 TO 1000000 DEFAULT 5 (Minimum number of genes required to be in a pathway.)
# PARAMETER conditional.testing: conditional.testing TYPE [yes: yes, no: no] DEFAULT yes (Conditional testing means that when a significant GO term is found, i.e. p-value is smaller than the specified thershold, that GO term is removed when testing the significance of its parent.)
# PARAMETER p.adjust.method: p.adjust.method TYPE [none: none, BH: BH, BY: BY] DEFAULT none (Method for adjusting the p-value in order to account for multiple testing. Because of the structure of GO, multiple testing is theoretically problematic, and using conditional.testing is a generally the preferred method. The correction can only be applied when no conditional.testing is performed.)
# PARAMETER over.or.under.representation: over.or.under.representation TYPE [over: over, under: under] DEFAULT over (Should over or under-represented classes be seeked?)

# Hypergeometrix test of gene enrichment to term categories
# Dario Greco 7.1.2007
# modified JTT 30.7.2007
# modified MG 5.2.2010
# modified IS 1.10.2010
# modified EK 2.7.2012

# Parameter settings (default) for testing purposes
#ontology<-"biological_process"
#p.value.threshold<-0.05
#minimum.population<-5
#conditional.testing<-"yes"
#p.adjust.method<-"none"
#over.or.under.representation<-"over"

# Loads the libraries
library(GOstats)
library(R2HTML)

# Loads the normalized data
dat<-read.table("normalized.tsv", header=TRUE, sep="\t", row.names=1)

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

#Account for the fact that alternative CDF annotation packages for some Affymetrix
#arrays are no longer supported and should be replaced with their Bioconductor counterpart 
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


# check for conditional testing and multiple testing correction
if (conditional.testing == 'no') {
  conditional <- FALSE
} else {
  if (p.adjust.method != 'none')
    stop('CHIPSTER-NOTE: Multiple testing correction can be applied only when performing unconditional testing. Please set conditional.testing to no, or p.adjust.method to none. Usually the preferred method is to use conditional testing.')
  conditional <- TRUE
}

# define the output variable
output <- data.frame(total=integer(0), expected=numeric(0), observed=integer(0), p.value=numeric(0), description=character(0), ontology=character(0))

if (ontology == 'biological_process' || ontology == 'all') {
  params <- new('GOHyperGParams', geneIds=myids, annotation=annotpkg, ontology='BP', pvalueCutoff=p.value.threshold, conditional=conditional, testDirection=over.or.under.representation)
  go <- hyperGTest(params)
  go.table <- summary(go, pvalue=2)
  if (nrow(go.table)>0) {
    go.table$Pvalue <- p.adjust(go.table$Pvalue, method=p.adjust.method)
    go.table <- go.table[go.table$Pvalue <= p.value.threshold & go.table$Size >= minimum.population,]
    if (nrow(go.table)>0) {
      rownames(go.table) <- go.table[,1]
      go.table <- go.table[,c(6, 4, 5, 2, 7)]
      go.table$ontology <- 'biological process'
      colnames(go.table) <- colnames(output)
      output <- rbind(output, go.table)
      go.table$description <- paste('<a href="http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=', rownames(go.table), '">', go.table$description, '</a>', sep='')
      HTML(go.table, file='hypergeo.html', append=TRUE, Border=0, innerBorder=1)
    }
  }
}

if (ontology == 'molecular_function' || ontology == 'all') {
  params <- new('GOHyperGParams', geneIds=myids, annotation=annotpkg, ontology='MF', pvalueCutoff=p.value.threshold, conditional=conditional, testDirection=over.or.under.representation)
  go <- hyperGTest(params)
  go.table <- summary(go, pvalue=2)
  if (nrow(go.table)>0) {
    go.table$Pvalue <- p.adjust(go.table$Pvalue, method=p.adjust.method)
    go.table <- go.table[go.table$Pvalue <= p.value.threshold & go.table$Size >= minimum.population,]
    if (nrow(go.table)>0) {
      rownames(go.table) <- go.table[,1]
      go.table <- go.table[,c(6, 4, 5, 2, 7)]
      go.table$ontology <- 'molecular function'
      colnames(go.table) <- colnames(output)
      output <- rbind(output, go.table)
      go.table$description <- paste('<a href="http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=', rownames(go.table), '">', go.table$description, '</a>', sep='')
      HTML(go.table, file='hypergeo.html', append=TRUE, Border=0, innerBorder=1)
    }
  }
}

if (ontology == 'cellular_component' || ontology == 'all') {
  params <- new('GOHyperGParams', geneIds=myids, annotation=annotpkg, ontology='CC', pvalueCutoff=p.value.threshold, conditional=conditional, testDirection=over.or.under.representation)
  go <- hyperGTest(params)
  go.table <- summary(go, pvalue=2)
  if (nrow(go.table)>0) {
    go.table$Pvalue <- p.adjust(go.table$Pvalue, method=p.adjust.method)
    go.table <- go.table[go.table$Pvalue <= p.value.threshold & go.table$Size >= minimum.population,]
    if (nrow(go.table)>0) {
      rownames(go.table) <- go.table[,1]
      go.table <- go.table[,c(6, 4, 5, 2, 7)]
      go.table$ontology <- 'cellular component'
      colnames(go.table) <- colnames(output)
      output <- rbind(output, go.table)
      go.table$description <- paste('<a href="http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=', rownames(go.table), '">', go.table$description, '</a>', sep='')
      HTML(go.table, file='hypergeo.html', append=TRUE, Border=0, innerBorder=1)
    }
  }
}

# write outputs
write.table(output, file='hypergeo.tsv', quote=FALSE, sep='\t')
if (nrow(output)==0)
  write('<html>\n\t<body>\n\t\tNo significant results found!</br />\n\t</body>\n</html>', file='hypergeo.html')

# EOF
