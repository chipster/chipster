# ANALYSIS "aCGH tools (beta testing)"/"GO enrichment for copy number aberrations" (Performs a statistical test for enrichment of GO terms in frequently aberrated genes. The input should be the output from the module Convert called aCGH data from probes to genes.)
# INPUT GENE_EXPRS gene-aberrations.tsv
# OUTPUT hypergeo-go.tsv, hypergeo-go.html
# PARAMETER aberrations [losses, gains, amplifications] DEFAULT losses (Whether to test enrichment of GO terms in frequently lost, gained or amplified genes.)
# PARAMETER frequency.threshold DECIMAL DEFAULT 0.3 (The minimum proportion of samples containing the particular type of aberration.)
# PARAMETER p.value.threshold DECIMAL DEFAULT 0.05 (P-value threshold.)

# pathways-acgh-hyperg-go.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-08-30

# load packages
library(org.Hs.eg.db)
library(GOstats)

# read input
dat <- read.table('gene-aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

# convert list of reference genes from Ensembl to Entrez IDs
ensembl.to.entrez <- as.list(org.Hs.egENSEMBL2EG)
reference.genes <- unique(unlist(ensembl.to.entrez[rownames(dat)]))

# check that we have something (i.e. that input file was in fact Ensembl IDs)
if (length(reference.genes)==0)
  stop('CHIPSTER-NOTE: The input file should contain a list of Ensembl Gene IDs. Usually as a result of running the module Convert called aCGH data from probes to genes.')

# detect the frequency column to use
if (aberrations == 'amplifications') {
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

# define output variables
output <- data.frame(p.value=numeric(0), odds.ratio=numeric(0), expected.count=numeric(0), count=integer(0), size=integer(0), term=character(0), ontology=character(0), direction=character(0))
ontology.labels <- c('BP'='biological process', 'MF'='molecular function', 'CC'='cellular component')
direction.labels <- c('over'='over-represented', 'under'='under-represented')

# run all six tests
ontologies <- c('BP', 'MF', 'CC')
testDirections <- c('over', 'under')
for (ontology in ontologies) {
  for (testDirection in testDirections) {
    params <- new('GOHyperGParams', geneIds=selected.genes, universeGeneIds=reference.genes, annotation='org.Hs.eg.db', ontology=ontology, pvalueCutoff=p.value.threshold, conditional=TRUE, testDirection=testDirection)
    go <- hyperGTest(params)
    go.table <- summary(go)
    if (nrow(go.table)>0) {
      rownames(go.table) <- go.table[,1]
      go.table[,1] <- NULL
      go.table$ontology <- ontology.labels[ontology]
      go.table$direction <- direction.labels[testDirection]
      colnames(go.table) <- colnames(output)
      output <- rbind(output, go.table)
      htmlReport(go, file='hypergeo-go.html', append=TRUE)
    }
  }
}

# write outputs
write.table(output, file='hypergeo-go.tsv', quote=FALSE, sep='\t')
if (nrow(output)==0)
  write('<html>\n\t<body>\n\t\tNo significant results found!</br />\n\t</body>\n</html>', file='hypergeo-go.html')

# EOF