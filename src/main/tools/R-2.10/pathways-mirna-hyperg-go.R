# ANALYSIS Pathways/"GO enrichment for miRNA targets" (Performs a statistical test for enrichments of GO terms in the predicted gene targets of a list of miRNA ID:s.)
# INPUT GENE_EXPRS normalized.tsv
# OUTPUT hyperg_go.tsv, hyperg_go.html
# PARAMETER ontology [all, biological_process, molecular_function, cellular_component] DEFAULT biological_process (The ontology to be analyzed.)
# PARAMETER p.value.threshold DECIMAL DEFAULT 0.05 (P-value threshold.)
# PARAMETER minimum.population INTEGER FROM 1 TO 1000000 DEFAULT 2 (Minimum number of genes required to be in a pathway.)
# PARAMETER conditional.testing [yes, no] (Conditional testing means that when a significant GO term is found, i.e. p-value is smaller than the specified thershold, that GO term is removed when testing the significance of its parent.)
# PARAMETER p.adjust.method [none, BH, BY] DEFAULT none (Method for adjusting the p-value in order to account for multiple testing. Because of the structure of GO, multiple testing is theoretically problematic, and using conditional.testing is a generally the preferred method. The correction can only be applied when no conditional.testing is performed.)
# PARAMETER over.or.under.representation [over, under] DEFAULT over (Should over or under-represented classes be seeked?)
# PARAMETER species [human, mouse, rat] DEFAULT human (The species for which the miRNA:s have been analyzed.)

# miRNA hypergeometric test for GO
# MG, 4.11.2009
# IS, 16.9.2010 rewritten to use intersection of pictar/targetscan for predicted targets, and GOstats for hypergeometric testing

# load packages
library(RmiR.Hs.miRNA) # what about other species ???
library(GOstats)
library(R2HTML)

# read input
dat <- read.table('normalized.tsv', header=TRUE, sep='\t', row.names=1)

# extracts identifiers
mirna_ids <- as.character(rownames(dat))

# load target predictions from pictar and targescan, intersect to build list of reference genes
pictar <- dbReadTable(RmiR.Hs.miRNA_dbconn(), 'pictar')[,1:2]
targetscan <- dbReadTable(RmiR.Hs.miRNA_dbconn(), 'targetscan')[,1:2]
reference.genes <- unique(intersect(pictar$gene_id, targetscan$gene_id))

# pick targets of the specified miRNAs
pictar <- pictar[pictar[,1] %in% mirna_ids,]
targetscan <- targetscan[targetscan[,1] %in% mirna_ids,]
selected.genes <- unique(intersect(pictar$gene_id, targetscan$gene_id))

# check for conditional testing and multiple testing correction
if (conditional.testing == 'no') {
	conditional <- FALSE
} else {
	if (p.adjust.method != 'none')
		stop('CHIPSTER-NOTE: Multiple testing correction can be applied only when performing unconditional testing. Please set conditional.testing to no, or p.adjust.method to none. Usually the preferred method is to use conditional testing.')
	conditional <- TRUE
}

# define annotation package
if (species == 'rat') {
	annotpkg <- 'org.Rn.eg.db'
} else if (species == 'mouse') {
	annotpkg <- 'org.Mm.eg.db'
} else {
	annotpkg <- 'org.Hs.eg.db'
}

# define the output variable
output <- data.frame(total=integer(0), expected=numeric(0), observad=integer(0), p.value=numeric(0), description=character(0), ontology=character(0))

if (ontology == 'biological_process' || ontology == 'all') {
	params <- new('GOHyperGParams', geneIds=selected.genes, universeGeneIds=reference.genes, annotation=annotpkg, ontology='BP', pvalueCutoff=p.value.threshold, conditional=conditional, testDirection=over.or.under.representation)
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
			HTML(go.table, file='hyperg_go.html', append=TRUE, Border=0, innerBorder=1)
		}
	}
}

if (ontology == 'molecular_function' || ontology == 'all') {
	params <- new('GOHyperGParams', geneIds=selected.genes, universeGeneIds=reference.genes, annotation=annotpkg, ontology='MF', pvalueCutoff=p.value.threshold, conditional=conditional, testDirection=over.or.under.representation)
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
			HTML(go.table, file='hyperg_go.html', append=TRUE, Border=0, innerBorder=1)
		}
	}
}

if (ontology == 'cellular_component' || ontology == 'all') {
	params <- new('GOHyperGParams', geneIds=selected.genes, universeGeneIds=reference.genes, annotation=annotpkg, ontology='CC', pvalueCutoff=p.value.threshold, conditional=conditional, testDirection=over.or.under.representation)
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
			HTML(go.table, file='hyperg_go.html', append=TRUE, Border=0, innerBorder=1)
		}
	}
}

# write outputs
write.table(output, file='hyperg_go.tsv', quote=FALSE, sep='\t')
if (nrow(output)==0)
	write('<html>\n\t<body>\n\t\tNo significant results found!</br />\n\t</body>\n</html>', file='hyperg_go.html')

# EOF