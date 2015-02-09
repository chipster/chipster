# TOOL pathways-ngs-hyperg-go.R: "GO enrichment for list of genes" (Performs a statistical test for enrichment of GO terms in a query list of genes. The input file should be from the tool \"Find unique and annotated genes\".)
# INPUT gene-list.tsv: gene-list.tsv TYPE GENERIC 
# OUTPUT hypergeo-go.tsv: hypergeo-go.tsv 
# OUTPUT hypergeo-go.html: hypergeo-go.html 
# PARAMETER genome: "Genome" TYPE [org.Hs.eg.db: "human", org.Mm.eg.db: "mouse", org.Rn.eg.db: "rat", org.Ag.eg.db: "arabidopsis", org.Dm.eg.db: "fruitfly"] DEFAULT org.Hs.eg.db (Reference organism.)
# PARAMETER column: "Identifier column" TYPE COLUMN_SEL DEFAULT entrezgene (Column containing gene identifiers. Identifiers should either be Entrez or EnsEMBL identifiers. When left blank, the first column is used.)
# PARAMETER OPTIONAL ontology: Ontology TYPE [all: all, biological_process: "biological process", molecular_function: "molecular function", cellular_component: "cellular component"] DEFAULT biological_process (The ontology to be analyzed.)
# PARAMETER OPTIONAL p.value.threshold: "P-value threshold" TYPE DECIMAL DEFAULT 0.05 (P-value threshold.)
# PARAMETER OPTIONAL minimum.population: "Minimum population" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 5 (Minimum number of genes required to be in a pathway.)
# PARAMETER OPTIONAL conditional.testing: "Conditional testing" TYPE [yes: yes, no: no] DEFAULT yes (Conditional testing means that when a significant GO term is found, i.e. p-value is smaller than the specified thershold, that GO term is removed when testing the significance of its parent.)
# PARAMETER OPTIONAL p.adjust.method: "P-value adjustment method" TYPE [none: none, BH: BH, BY: BY] DEFAULT none (Method for adjusting the p-value in order to account for multiple testing. Because of the structure of GO, multiple testing is theoretically problematic, and using conditional.testing is a generally the preferred method. The correction can only be applied when no conditional.testing is performed.)
# PARAMETER OPTIONAL over.or.under.representation: "Over or under-representation" TYPE [over: over, under: under] DEFAULT over (Should over or under-represented classes be seeked?)

# pathways-ngs-hyperg-go.R
# 12.02.2010 MG, Created
# 09.05.2014 MK, Added several new organisms. Added possibility to annotate Ensembl IDs as well.
# 25.06.2014 EK, Clarified parameters.

# load packages
library(genome, character.only=T)
annotpkg <- gsub(".db", "", genome)
library(GOstats)
library(R2HTML)

# read input
dat <- read.table('gene-list.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

# convert list of reference genes from Ensembl to Entrez IDs
# ensembl.to.entrez <- as.list(org.Hs.egENSEMBL2EG)
lib <- paste(annotpkg, "ENSEMBL2EG", sep="")
ensembl.to.entrez <- as.list(get(lib))
reference.genes <- unique(unlist(ensembl.to.entrez))
ens_identifiers <- names(unlist(ensembl.to.entrez))

# see if column is the first "empty" element, indicating that the user wishes to use row.names
if(length(column) == 0 || column == " " || column == "") {
	selected.genes <- row.names(dat)
} else {
	column <- grep(paste("^", column, "$", sep=""), colnames(dat))
	selected.genes <- as.character(dat[,column])
}

# check if IDs match more EnsEMBL or Entrez genes
if(length(intersect(selected.genes, reference.genes)) < length(intersect(selected.genes, ens_identifiers))) {
	selected.genes <- unique(unlist(ensembl.to.entrez[selected.genes]))
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
	params <- new('GOHyperGParams', geneIds=selected.genes, universeGeneIds=reference.genes, annotation=genome, ontology='BP', pvalueCutoff=p.value.threshold, conditional=conditional, testDirection=over.or.under.representation)
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
	params <- new('GOHyperGParams', geneIds=selected.genes, universeGeneIds=reference.genes, annotation=genome, ontology='MF', pvalueCutoff=p.value.threshold, conditional=conditional, testDirection=over.or.under.representation)
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
	params <- new('GOHyperGParams', geneIds=selected.genes, universeGeneIds=reference.genes, annotation=genome, ontology='CC', pvalueCutoff=p.value.threshold, conditional=conditional, testDirection=over.or.under.representation)
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
