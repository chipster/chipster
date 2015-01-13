# TOOL pathways-mirna-hyperg-go.R: "GO enrichment for miRNA targets" (Given a list of miRNA identifiers, tests for enrichment of GO terms in their predicted gene targets.)
# INPUT normalized.tsv: normalized.tsv TYPE GENERIC 
# OUTPUT hyperg_go.tsv: hyperg_go.tsv 
# OUTPUT hyperg_go.html: hyperg_go.html
# PARAMETER species: Organism TYPE [human: human, mouse: mouse, rat: rat] DEFAULT human (From which organism does the data come from)
# PARAMETER OPTIONAL ontology: "GO ontology to be tested" TYPE [all: all, biological_process: "biological process", molecular_function: "molecular function", cellular_component: "cellular component"] DEFAULT biological_process (The ontology to be analyzed.)
# PARAMETER OPTIONAL p.value.threshold: "p-value threshold" TYPE DECIMAL DEFAULT 0.05 (P-value threshold.)
# PARAMETER OPTIONAL minimum.population: "Minimum population" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 5 (Minimum number of genes required to be in a pathway.)
# PARAMETER OPTIONAL conditional.testing: "Use conditional testing" TYPE [yes: yes, no: no] DEFAULT yes (Conditional testing means that when a significant GO term is found, i.e. p-value is smaller than the specified threshold, that GO term is removed when testing the significance of its parent.)
# PARAMETER OPTIONAL p.adjust.method: "p-value adjustment method" TYPE [none: none, BH: BH, BY: BY] DEFAULT none (Method for adjusting the p-value in order to account for multiple testing. Because of the structure of GO, multiple testing is theoretically problematic, and using conditional.testing is a generally the preferred method. The correction can only be applied when no conditional.testing is performed.)
# PARAMETER OPTIONAL over.or.under.representation: "Over- or underrepresentation" TYPE [over: over, under: under] DEFAULT over (Should over or under-represented classes be seeked?)
# PARAMETER OPTIONAL database: "Database for human miRNA targets" TYPE [PicTar: PicTar, TargetScan: TargetScan, both: both] DEFAULT both (For human data this parameter defines whether to fetch predicted gene targets from either PicTar or TargetScan database, or whether to restrict the enrichment analysis to those targets that are common to both databases. The last option is recommended to minimize false positives.)

# miRNA hypergeometric test for GO
# MG, 4.11.2009
# IS, 1.10.2010, rewritten to use GOstats
# MG, 7.3.2011, added the "database" parameter
# EK and MK, 10.9.2013, miRNA names are forced to lower case so that they work with the local miRNA_mappings files (used for mouse and rat)
# MK, remember to modify pathways-mirna-hyperg-kegg. Mm uses now targetscan predictions

# load packages
library(GOstats)
library(R2HTML)

# read input
dat <- read.table('normalized.tsv', header=TRUE, sep='\t', row.names=1)

# extracts identifiers
if("id" %in% colnames(dat)) {
	mirna_ids <- as.character(dat$id)
} else {
	mirna_ids <- as.character(rownames(dat))
}

# If convert genomic BAM file has been used, table has a column which name is sequence
if("sequence" %in% colnames(dat) && (length(grep(dat$sequence[1], mirna_ids[1])) > 0)) {
	mirna_id_list <- strsplit(as.vector(mirna_ids), "_")
	mirna_ids <- NULL
	for(i in 1:length(mirna_id_list)) {
		#remove last three sections of the id
		mirna_ids <- c(mirna_ids, paste(unlist(mirna_id_list[i])[1:((length(unlist(mirna_id_list[i])))-3)], collapse="_"))
	}
}

mirna_ids <- tolower(mirna_ids)

# check for conditional testing and multiple testing correction
if (conditional.testing == 'no') {
	conditional <- FALSE
} else {
	if (p.adjust.method != 'none')
		stop('CHIPSTER-NOTE: Multiple testing correction can be applied only when performing unconditional testing. Please set conditional.testing to no, or p.adjust.method to none. Usually the preferred method is to use conditional testing.')
	conditional <- TRUE
}

if (species == 'mouse') {
	if (database == "PicTar" || database == "both") {
		stop('CHIPSTER-NOTE: Only TargetScan database is available for mouse')
	}

	library(targetscan.Mm.eg.db)
	
	reference.genes <- unique(names(as.list(targetscan.Mm.egTARGETS)))
	selected.genes <- NULL
	for(i in 1:length(mirna_ids)) {
		mirna.fam <- try(mget(as.vector(mirna_ids[i]), revmap(targetscan.Mm.egFAMILY2MIRBASE)), silent=T)
		if(class(mirna.fam) != "try-error") {
			mirna.targets <- try(mget(unlist(mirna.fam), revmap(targetscan.Mm.egTARGETS)), silent=T)
			if(class(mirna.targets) != "try-error") {	
				selected.genes <- c(selected.genes, unique(unlist(mirna.targets)))
			}
		}
	}

	#for(i in 1:length(mirna_ids)) {
	#	a <- try(mget(mirna_ids[i], revmap(targetscan.Mm.egTARGETS)), silent=T)
	#	if(class(a) != "try-error") {
	#		selected.genes <- c(selected.genes, unlist(a))
	#	}
	#}
	#selected.genes <- unique(selected.genes)

	if (length (selected.genes) == 0) {
		stop("CHIPSTER-NOTE: No target genes were found for the input list of miRNA names. Please make sure that you are using official miRNA names.")
	}  
	annotpkg <- 'org.Mm.eg.db'
	
} else if (species == 'rat') {
	library(org.Rn.eg.db)
	# targets <- miRBase2df.fun(url="ftp://ftp.sanger.ac.uk/pub/mirbase/targets/v5/arch.v5.txt.rattus_norvegicus.zip")
	targets <- read.table(file.path(chipster.tools.path, "miRNA_mappings", "mirna_mappings_rnorvegicus.txt"), sep="\t")
	
	ensembl.to.entrez <- as.list(org.Rn.egENSEMBLTRANS2EG)
	reference.genes <- unique(unlist(ensembl.to.entrez[unique(targets$tran)]))
	selected.genes <- unique(unlist(ensembl.to.entrez[unique(targets[tolower(targets$mir) %in% mirna_ids, 'tran'])]))
	
	# check that it was indeed possible to identify targets for the input list of miRNA names
	if (length (selected.genes) == 0) {
		stop("CHIPSTER-NOTE: No target genes were found for the input list of miRNA names. Please make sure that you are using official miRNA names.")
	}  
	annotpkg <- 'org.Rn.eg.db'
	
} else {
	library(RmiR.Hs.miRNA)
	library(org.Hs.eg.db)
	# load target predictions from pictar and targescan, intersect to build list of reference genes

	reference.genes <- NULL
	selected.genes <- NULL
	if (database == "PicTar" || database == "both") {
		pictar <- dbReadTable(RmiR.Hs.miRNA_dbconn(), 'pictar')[,1:2]
		reference.genes <- c(reference.genes, unique(pictar$gene_id))

		pictar <- pictar[tolower(pictar[,1]) %in% mirna_ids,]
		selected.genes <- c(selected.genes, unique(pictar$gene_id))
	}
	if (database == "TargetScan" || database == "both") {
		targetscan <- dbReadTable(RmiR.Hs.miRNA_dbconn(), 'targetscan')[,1:2]
		reference.genes <- c(reference.genes, unique(targetscan$gene_id))

		targetscan <- targetscan[tolower(targetscan[,1]) %in% mirna_ids,]
		selected.genes <- c(selected.genes, unique(targetscan$gene_id))
	}
	
	# check that it was indeed possible to identify targets for the input list of miRNA names
	if (length (selected.genes) == 0) {
		stop("CHIPSTER-NOTE: No target genes were found for the input list of miRNA names. Please make sure that you are using official miRNA names.")
	}  
	annotpkg <- 'org.Hs.eg.db'
}

# define the output variable
output <- data.frame(total=integer(0), expected=numeric(0), observed=integer(0), p.value=numeric(0), description=character(0), ontology=character(0))

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