# TOOL annotate-ensembl-ids.R: "Annotate Ensembl-IDs" (Annotates Ensembl IDs with gene symbols and descriptions, and adds the results to the datafile.)
# INPUT genelist.tsv: genelist.tsv TYPE GENERIC
# OUTPUT annotated.tsv: annotated.tsv 
# PARAMETER species: Species TYPE [human: human, mouse: mouse, rat: rat] DEFAULT human (The species needs to be specified in order to map genes to the genomic coordinates.)


# ML, 05.02.2015

# Toimii, mut vain human, mouse, rat. / Testattu vasta human ja htseq-count.R:n tulostiedosto. 
# Tarvitaan vielä parametri ja siitä riippuva IF että onko eka sarake vai rivinimet ensemblID.
# Kestää aika kauan.

# Loads libraries into memory
library(biomaRt)

# Loads the data
file <- c("genelist.tsv")
#dat <- read.table(file, header=T, sep="\t", row.names=1)
dat <- read.table(file, header=T, sep="\t")
genes <- dat$id
		# genes <- dat$id ???? Onko kaikissa nimetty näin?


# Fetch the gene symbols and descriptions from ENSEMBL using biomaRt
if (species=="human") {
	dataset <- "hsapiens_gene_ensembl"
	filt <- "hgnc_symbol"
}
if (species=="mouse") {
	dataset <- "mmusculus_gene_ensembl"
	filt <- "mgi_symbol"
}
if (species=="rat") {
	dataset <- "rnorvegicus_gene_ensembl"
	filt <- "rgd_symbol"
}

ensembl <- useMart("ensembl", dataset=dataset)
# annotated_genes <- getBM(mart=ensembl, attributes=c(filt,"chromosome_name","start_position","end_position"), filters=filt, values=gene_symbols)
# modified this:
# genes_ensembl_org <- getBM(mart=ensembl, attributes=c("entrezgene", "ensembl_gene_id", "external_gene_name", "description"), filters="ensembl_gene_id", values=genes)
genes_ensembl_org <- getBM(attributes <- c("entrezgene", "ensembl_gene_id", "external_gene_name", "description"), filters = "ensembl_gene_id", values = genes, mart = ensembl, uniqueRows=T)


#ensembl = useMart("ensembl")
#ensembl_hs = useDataset("hsapiens_gene_ensembl", mart=ensembl)


# genes <- dat$id ????
# genes_ensembl_org <- getBM(attributes <- c("entrezgene", "ensembl_gene_id", "external_gene_name", "description", "start_position", "end_position", "chromosome_name"), filters = "entrezgene", values = genes, mart = ensembl_hs, uniqueRows=T)
#genes_ensembl_org <- getBM(attributes <- c("entrezgene", "ensembl_gene_id", "external_gene_name", "description"), filters = "ensembl_gene_id", values = genes, mart = ensembl_hs, uniqueRows=T)

pmatch_table		<- pmatch(genes, genes_ensembl_org[,2], duplicates.ok=T)

ensembl_table		<- as.data.frame(matrix(NA, nrow=length(genes), ncol=8))
ensembl_table[which(!is.na(pmatch_table)),] <- genes_ensembl_org[pmatch_table[(!is.na(pmatch_table))], ];
rownames(ensembl_table)	<- genes;
colnames(ensembl_table) <- colnames(genes_ensembl_org);

results <- cbind(ensembl_table[,3:4], dat);
names(results) <- c("External gene name", "Description", names(dat))

# write result table to output
#results <- cbind(dat, genes);
write.table(results, file="annotated.tsv", col.names=T, quote=F, sep="\t", row.names=F)
#write.table(dat, file="htseq-counts.tsv", col.names=T, quote=F, sep="\t", row.names=F)



