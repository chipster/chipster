# ANALYSIS "aCGH tools (beta testing)"/"GO enrichment for copy number aberrations" (Performs a statistical test for enrichment of GO terms in frequently aberrated genes. The input should be the output from the module Convert called aCGH data from probes to genes.)
# INPUT GENE_EXPRS gene-aberrations.tsv
# OUTPUT hypergeo-go.tsv, hypergeo-go.html
# PARAMETER aberrations [losses, gains, amplifications, gains_and_amplifications, any_aberration] DEFAULT losses (Whether to test enrichment of GO terms in frequently lost, gained or amplified genes.)
# PARAMETER frequency.threshold DECIMAL DEFAULT 0.5 (The minimum proportion of samples containing the particular type of aberration.)
# PARAMETER ontology [all, biological_process, molecular_function, cellular_component] DEFAULT biological_process (The ontology to be analyzed.)
# PARAMETER p.value.threshold DECIMAL DEFAULT 0.05 (P-value threshold.)
# PARAMETER minimum.population INTEGER FROM 1 TO 1000000 DEFAULT 2 (Minimum number of genes required to be in a pathway.)
# PARAMETER conditional.testing [yes, no] (Conditional testing means that when a significant GO term is found, i.e. p-value is smaller than the specified thershold, that GO term is removed when testing the significance of its parent.)
# PARAMETER p.adjust.method [none, BH, BY] DEFAULT none (Method for adjusting the p-value in order to account for multiple testing. Because of the structure of GO, multiple testing is theoretically problematic, and using conditional.testing is a generally the preferred method. The correction can only be applied when no conditional.testing is performed.)
# PARAMETER over.or.under.representation [over, under] DEFAULT over (Should over or under-represented classes be seeked?)

# pathways-acgh-hyperg-go.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-09-16

# load packages
library(org.Hs.eg.db)
library(GOstats)
library(R2HTML)

# read input
dat <- read.table('gene-aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

# convert list of reference genes from Ensembl to Entrez IDs
ensembl.to.entrez <- as.list(org.Hs.egENSEMBL2EG)
reference.genes <- unique(unlist(ensembl.to.entrez[rownames(dat)]))

# check that we have something (i.e. that input file was in fact Ensembl IDs)
if (length(reference.genes)==0)
  stop('CHIPSTER-NOTE: The input file should contain a list of Ensembl Gene IDs. Usually as a result of running the module Convert called aCGH data from probes to genes.')

# detect the frequency column to use
if (aberrations == 'any_aberration') {
  if ('amp.freq' %in% colnames(dat))
    dat$gain.freq <- dat$gain.freq + dat$amp.freq
  dat$loss.freq <- dat$loss.freq + dat$gain.freq
  column <- 'loss.freq'
} else if (aberrations == 'gains_and_amplifications') {
  if ('amp.freq' %in% colnames(dat))
    dat$gain.freq <- dat$gain.freq + dat$amp.freq
  column <- 'gain.freq'
} else if (aberrations == 'amplifications') {
  column <- 'amp.freq'
} else if (aberrations == 'gains') {
  column <- 'gain.freq'
} else {
  column <- 'loss.freq'
}

# check that the frequency column exists
if (!column %in% colnames(dat))
  stop('CHIPSTER-NOTE: The required frequency column not found in file: ', column)

# extract the list of frequently aberrated genes
selected.genes <- unique(unlist(ensembl.to.entrez[rownames(dat[dat[,column] >= frequency.threshold,])]))

# check that we have a list of genes to test
if (length(reference.genes)==0)
  stop('CHIPSTER-NOTE: There were no aberrated genes above the selected threshold (', frequency.threshold, '). Please choose a lower threshold.')

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
  params <- new('GOHyperGParams', geneIds=selected.genes, universeGeneIds=reference.genes, annotation='org.Hs.eg.db', ontology='BP', pvalueCutoff=p.value.threshold, conditional=conditional, testDirection=over.or.under.representation)
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
      HTML(go.table, file='hypergeo-go.html', append=TRUE, Border=0, innerBorder=1)
    }
  }
}

if (ontology == 'molecular_function' || ontology == 'all') {
  params <- new('GOHyperGParams', geneIds=selected.genes, universeGeneIds=reference.genes, annotation='org.Hs.eg.db', ontology='MF', pvalueCutoff=p.value.threshold, conditional=conditional, testDirection=over.or.under.representation)
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
      HTML(go.table, file='hypergeo-go.html', append=TRUE, Border=0, innerBorder=1)
    }
  }
}

if (ontology == 'cellular_component' || ontology == 'all') {
  params <- new('GOHyperGParams', geneIds=selected.genes, universeGeneIds=reference.genes, annotation='org.Hs.eg.db', ontology='CC', pvalueCutoff=p.value.threshold, conditional=conditional, testDirection=over.or.under.representation)
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
      HTML(go.table, file='hypergeo-go.html', append=TRUE, Border=0, innerBorder=1)
    }
  }
}

# write outputs
write.table(output, file='hypergeo-go.tsv', quote=FALSE, sep='\t')
if (nrow(output)==0)
  write('<html>\n\t<body>\n\t\tNo significant results found!</br />\n\t</body>\n</html>', file='hypergeo-go.html')

# EOF