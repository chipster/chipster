# TOOL acgh-pathways-go.R: "GO enrichment for called gene copy numbers" (Performs a statistical test for enrichment of GO terms in frequently aberrated genes. The input should be the output from the tool Detect genes from called copy number data.)
# INPUT gene-aberrations.tsv: gene-aberrations.tsv TYPE GENE_EXPRS 
# OUTPUT hypergeo-go.tsv: hypergeo-go.tsv 
# OUTPUT hypergeo-go.html: hypergeo-go.html 
# PARAMETER aberrations: Aberrations TYPE [all_aberrations: "all aberrations", losses_and_deletions: "losses and deletions", deletions: deletions, losses: losses, gains: gains, amplifications: amplifications, gains_and_amplifications: "gains and amplifications"] DEFAULT all_aberrations (Whether to test enrichment of GO terms in frequently lost, gained or amplified genes.)
# PARAMETER frequency.threshold: "Frequency threshold" TYPE DECIMAL DEFAULT 0.5 (The minimum proportion of samples containing the particular type of aberration.)
# PARAMETER ontology: Ontology TYPE [all: all, biological_process: "biological process", molecular_function: "molecular function", cellular_component: "cellular component"] DEFAULT biological_process (The ontology to be analyzed.)
# PARAMETER p.value.threshold: "p-value threshold" TYPE DECIMAL DEFAULT 0.05 (P-value threshold.)
# PARAMETER minimum.population: "Minimum population" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 2 (Minimum number of genes required to be in a pathway.)
# PARAMETER conditional.testing: "Conditional testing" TYPE [yes: yes, no: no] DEFAULT yes (Conditional testing means that when a significant GO term is found, i.e. p-value is smaller than the specified thershold, that GO term is removed when testing the significance of its parent.)
# PARAMETER p.adjust.method: "p-adjust method" TYPE [none: none, BH: BH, BY: BY] DEFAULT none (Method for adjusting the p-value in order to account for multiple testing. Because of the structure of GO, multiple testing is theoretically problematic, and using conditional.testing is a generally the preferred method. The correction can only be applied when no conditional.testing is performed.)
# PARAMETER over.or.under.representation: "Over- or under representation" TYPE [over: over, under: under] DEFAULT over (Should over or under-represented classes be seeked?)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-24

source(file.path(chipster.common.path, 'library-Chipster.R'))
library(org.Hs.eg.db)
library(GOstats)
library(R2HTML)

# read input
dat <- readData("gene-aberrations.tsv")

# detect the frequency column to use
if (aberrations == 'gains_and_amplifications') {
  if ('amp.freq' %in% colnames(dat))
    dat$gain.freq <- dat$gain.freq + dat$amp.freq
  column <- 'gain.freq'
} else if (aberrations == 'amplifications') {
  column <- 'amp.freq'
} else if (aberrations == 'gains') {
  column <- 'gain.freq'
} else if (aberrations == 'losses'){
  column <- 'loss.freq'
} else if (aberrations == 'deletions'){
  column <- 'del.freq'
} else if (aberrations == 'losses_and_deletions') {
  if ('del.freq' %in% colnames(dat))
    dat$loss.freq <- dat$loss.freq + dat$del.freq
  column <- 'loss.freq'
} else {
  if ('del.freq' %in% colnames(dat))
    dat$loss.freq <- dat$loss.freq + dat$del.freq
  if ('amp.freq' %in% colnames(dat))
    dat$gain.freq <- dat$gain.freq + dat$amp.freq
  dat$loss.freq <- dat$loss.freq + dat$gain.freq
  column <- 'loss.freq'
}

# check that the frequency column exists
if (!column %in% colnames(dat))
  stop('CHIPSTER-NOTE: The required frequency column not found in file: ', column)

if ("entrez" %in% colnames(dat)) {
  # new version of acgh-convert-from-probes-to-genes.R contains column "entrez"
  reference.genes <- dat$entrez
  selected.genes <- dat[dat[,column] >= frequency.threshold, "entrez"]
} else {
  # old version of acgh-convert-from-probes-to-genes.R contained ensembl ids
  ensembl.to.entrez <- as.list(org.Hs.egENSEMBL2EG)
  reference.genes <- unique(unlist(ensembl.to.entrez[rownames(dat)]))
  selected.genes <- unique(unlist(ensembl.to.entrez[rownames(dat[dat[,column] >= frequency.threshold,])]))
}

if (length(reference.genes) == 0)
  stop('CHIPSTER-NOTE: No gene IDs found in the input file, please first run the tool Detect genes from called copy number data.')

# check that we have a list of genes to test
if (length(selected.genes) == 0)
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
  if (nrow(go.table) > 0) {
    go.table$Pvalue <- p.adjust(go.table$Pvalue, method=p.adjust.method)
    go.table <- go.table[go.table$Pvalue <= p.value.threshold & go.table$Size >= minimum.population,]
    if (nrow(go.table) > 0) {
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
  if (nrow(go.table) > 0) {
    go.table$Pvalue <- p.adjust(go.table$Pvalue, method=p.adjust.method)
    go.table <- go.table[go.table$Pvalue <= p.value.threshold & go.table$Size >= minimum.population,]
    if (nrow(go.table) > 0) {
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
  if (nrow(go.table) > 0) {
    go.table$Pvalue <- p.adjust(go.table$Pvalue, method=p.adjust.method)
    go.table <- go.table[go.table$Pvalue <= p.value.threshold & go.table$Size >= minimum.population,]
    if (nrow(go.table) > 0) {
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

output$expected <- signif(output$expected, digits=3)
output$p.value <- signif(output$p.value, digits=3)

# write outputs
writeData(output, "hypergeo-go.tsv")
if (nrow(output) == 0)
  write('<html>\n\t<body>\n\t\tNo significant results found!</br />\n\t</body>\n</html>', file='hypergeo-go.html')

# EOF
