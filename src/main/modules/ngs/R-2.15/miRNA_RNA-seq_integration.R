# TOOL miRNA_RNA-seq_integration.R: "Correlate miRNA-seq and RNA-seq data" (Detects miRNA target genes whose expression correlates with miRNA expression, either negatively or positively. Note that you need miRNA-seq and RNA-seq data from the same samples. The matching pairs need to be indicated with numbers in phenodata. This tool works only for human data currently.)
# INPUT normalized_mirna.tsv: "miRNA expression table" TYPE GENE_EXPRS 
# INPUT normalized_gene.tsv: "RNA expression table" TYPE GENE_EXPRS 
# INPUT phenodata_mirna.tsv: "Phenodata for miRNA" TYPE GENERIC 
# INPUT phenodata_gene.tsv: "Phenodata for RNA" TYPE GENERIC 
# OUTPUT OPTIONAL full_correlation_matrix.tsv: full_correlation_matrix.tsv
# OUTPUT OPTIONAL correlation_annotated_and_expressed_miRNAs.tsv: correlation_annotated_and_expressed_miRNAs.tsv
# OUTPUT OPTIONAL correlation_known_interactions_only.tsv: correlation_known_interactions_only.tsv
# PARAMETER order.column.mirna: "Phenodata column indicating sample order in miRNA data" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing the order of the samples with numbers, so that RNA and miRNA data can be correctly matched.)
# PARAMETER order.column.gene: "Phenodata column indicating sample order in RNA data" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing the order of the samples with numbers, so that RNA and miRNA data can be correctly matched.)
# PARAMETER OPTIONAL normalization.method: "Normalization" TYPE [none, cpm, TMM] DEFAULT none (Should the miRNA and RNA counts be normalized. Available methods are counts per million (cpm\) and trimmed mean of M-values (TMM\) based on the edgeR package.)
# PARAMETER OPTIONAL filtering.method: "Filter by" TYPE [correlation, p.value] DEFAULT correlation (Should miRNA-RNA pairs be filtered by correlation or by p-value.)
# PARAMETER OPTIONAL filter.threshold: "Filtering threshold" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.90 (Filtering cut-off.)
# PARAMETER OPTIONAL save.full.matrix: "Output also the full miRNA-RNA correlation matrix" TYPE [yes, no] DEFAULT no (This (large\) matrix contains correlations between all miRNAs and genes, with no filtering applied.)


# Correlation analysis of miRNA-seq and RNA-seq data
# JTT 2013-08-24

## setwd("C:\\Users\\Jarno Tuimala\\Desktop\\Chipster2013\\miRNA_rna-seq\\sample_data")
## data_1<-read.table(file="normalized-mirna.tsv", header=T, sep="\t", row.names=1)
## data_2<-read.table(file="normalized-mrna.tsv", header=T, sep="\t", row.names=1)
## phenodata_1 <- read.table("phenodata-mirna.tsv", header=T, sep="\t")
## phenodata_2 <- read.table("phenodata-mrna.tsv", header=T, sep="\t")
## order.column.mirna<-"sample"
## order.column.gene<-"sample"
## filtering.method<-"correlation"
## filter.threshold<-0.90
## save.full.matrix<-"no"

# Loads the normalized data and phenodata files
data_1 <- read.table(file="normalized_mirna.tsv", header=T, sep="\t", row.names=1)
data_2 <- read.table(file="normalized_gene.tsv", header=T, sep="\t", row.names=1)
phenodata_1 <- read.table("phenodata_mirna.tsv", header=T, sep="\t")
phenodata_2 <- read.table("phenodata_gene.tsv", header=T, sep="\t")

# Figure out which is the miRNA data
if (phenodata_1$chiptype[1] == "miRNA") {
	mirna.phenodata <- phenodata_1
	mirna.data <- data_1
	gene.phenodata <- phenodata_2
	gene.data <- data_2
}
if (phenodata_2$chiptype[1] == "miRNA") {
	mirna.phenodata <- phenodata_2
	mirna.data <- data_2
	gene.phenodata <- phenodata_1
	gene.data <- data_1
}

# Separates expression values and other columns
mirna.data.2 <- mirna.data[,grep("chip", names(mirna.data))]
gene.data.2 <- gene.data[,grep("chip", names(gene.data))]

# pick those samples that do have a matching pair
common.samples <- intersect(mirna.phenodata[,order.column.mirna], gene.phenodata[,order.column.gene])
rownames(mirna.phenodata) <- mirna.phenodata[,order.column.mirna]
rownames(gene.phenodata) <- gene.phenodata[,order.column.gene]
mirna.phenodata$n <- 1:nrow(mirna.phenodata)
gene.phenodata$n <- 1:nrow(gene.phenodata)
mirna.order <- mirna.phenodata[common.samples, 'n']
gene.order <- gene.phenodata[common.samples, 'n']

# Arrange the columns in the two data sets so that they match
mirna.data.3 <- mirna.data.2[,order(mirna.order)]
gene.data.3 <- gene.data.2[,order(gene.order)]

# Normalization
if(normalization.method=="tmm") {
   library(edgeR)
   mirna3<-DGEList(mirna.data.3)
   gene3<-DGEList(gene.data.3)
   mirna3.1<-calcNormFactors(mirna3)
   gene3.1<-calcNormFactors(gene3)
   mirna3.2 = estimateCommonDisp(mirna3.1, verbose=TRUE)
   gene3.2 = estimateCommonDisp(gene3.1, verbose=TRUE)
   mirna.data.3<-mirna3.2$pseudo.counts
   gene.data.3<-gene3.2$pseudo.counts
}
if(normalization.method=="cpm") {
   library(edgeR)
   mirna3<-DGEList(mirna.data.3)
   gene3<-DGEList(gene.data.3)
   mirna3.1<-calcNormFactors(mirna3)
   gene3.1<-calcNormFactors(gene3)
   mirna.data.3<-cpm(mirna3.1, normalized.lib.sizes=TRUE)
   gene.data.3<-cpm(gene3.1, normalized.lib.sizes=TRUE)
}

# Pearson correlation coefficients and the corresponding p-values are calculated for all possible miRNA-mRNA pairs
library(WGCNA)
corp<-function (x, y = NULL, use = "pairwise.complete.obs", alternative = c("two.sided", "less", "greater"), ...) {
    cor = cor(x, y, use = use, ...)
    x = as.matrix(x)
    finMat = !is.na(x)
    y = as.matrix(y)
    np = t(finMat) %*% (!is.na(y))
    T = sqrt(np - 2) * abs(cor)/sqrt(1 - cor^2)
    p = 2 * pt(T, np - 2, lower.tail = FALSE)
    list(cor = cor, p = p, nObs = np)
}
d<-corp(t(gene.data.3), t(mirna.data.3), use="pairwise.complete.obs")

# Convert possible ENSEMBL IDs to Entrez Gene
if(length(grep("ENS", rownames(gene.data)))>0) {
   id<-as.character(rownames(gene.data))
   library(org.Hs.eg.db)
   xx <- as.list(org.Hs.egENSEMBL2EG)
   dd<-as.data.frame(unlist(xx))
   id2<-as.data.frame(id)
   m<-merge(id2, dd, by.x="id", by.y="row.names", sort=F, all.x=T)
} else {
   m<-data.frame(id=rownames(gene.data), entrez.gene=rownames(gene.data))
}
colnames(m)<-c("id", "gene")

# Get targets for miRNAs as Entrez gene ids
library(RmiR.Hs.miRNA)
miranda <- dbReadTable(RmiR.Hs.miRNA_dbconn(), "miranda")[,c(2,1)]
mirbase <- dbReadTable(RmiR.Hs.miRNA_dbconn(), "mirbase")[,1:2]
targetscan <- dbReadTable(RmiR.Hs.miRNA_dbconn(), "targetscan")[,1:2]
pictar <- dbReadTable(RmiR.Hs.miRNA_dbconn(), "pictar")[,1:2]
tarbase <- dbReadTable(RmiR.Hs.miRNA_dbconn(), "tarbase")[,1:2]
mid<-rbind(miranda, mirbase, targetscan, pictar, tarbase)
mid[,1]<-tolower(mid[,1])
mid2<-mid[!duplicated(mid),]

# Some cleaning
rm(miranda, mirbase, targetscan, pictar, tarbase)
gc()

# Keep only miRNAs that are expressed in at least one sample and have target gene annotation
mirna.ind<-unique(which(as.character(rownames(mirna.data)) %in% as.character(mid$mature_miRNA)))
gene.ind<-which(as.character(m$gene) %in% unique(mid$gene_id[which(as.character(mid$mature_miRNA) %in% as.character(rownames(mirna.data.3)))]))
if(length(mirna.ind)>=1 & length(gene.ind)>=1) {
   df<-list(cor=d$cor[gene.ind,mirna.ind,drop=F], p=d$p[gene.ind,mirna.ind,drop=F], nObs=d$nObs[gene.ind,mirna.ind,drop=F])
} else {
   stop("There are either no annotated miRNAs/mRNAs or any expressed miRNAs in your dataset! Aborting computations.")
}
col.ind<-colSums(df$cor, na.rm=TRUE)!=0
pval<-df$p[,col.ind,drop=F]
cval<-df$cor[,col.ind,drop=F]

# Write out results
if(save.full.matrix=="yes") {
   write.table(d, "full_correlation_matrix.tsv", col.names=T, row.names=T, sep="\t", quote=F)
}

# Some cleaning
rm(d, gene.data.2, mirna.data.2)
gc()

# Process the results
# Initiate the result table
ptemp<-c()
ctemp<-c()
for(i in 1:ncol(pval)) {
   ptemp<-c(ptemp, pval[,i])
   ctemp<-c(ctemp, cval[,i])
}

# Fill in the result table
res<-as.data.frame(matrix(ncol=6, nrow=length(ptemp), data=NA))
colnames(res)<-c("miRNA","entrez.gene.id","symbol","pearson.correlation.coefficient","p.value","original.id")
res[,1]<-rep(colnames(pval), each=nrow(pval))
res[,2]<-m$gene[match(names(ptemp), m$id)]
xx <- as.list(org.Hs.egSYMBOL)
dd<-as.data.frame(unlist(xx))
res[,3]<-dd[match(res[,2], rownames(dd)),1]
res[,4]<-as.vector(ctemp)
res[,5]<-as.vector(ptemp)
res[,6]<-names(ptemp)

# Some cleaning
rm(ptemp, ctemp)
gc()

# Filter the results on correlation or p-value
if(filtering.method=="correlation") {
   res2<-na.omit(res[res$pearson.correlation.coefficient>=filter.threshold | res$pearson.correlation.coefficient<=-filter.threshold,])
}

if(filtering.method=="p-value") {
   res2<-na.omit(res[res$p.value<=filter.threshold,])
}

# if(nrow(res2)==0) {
#    stop("No results left after threshold filtering! Aborting...")
# }

# Filter the result pairs on known interactions
mirnas<-unique(res2$miRNA)
res3<-c()
for(i in 1:length(mirnas)) {
	rt<-res2[res2$miRNA==mirnas[i],]
	mt<-mid2[mid2$mature_miRNA==mirnas[i],]
	res3<-rbind(res3,rt[rt$entrez.gene.id %in% mt$gene_id,])
}

#if(nrow(res3)==0) {
#   stop("No results left after known interactions filtering! Aborting...")
#}

# Write out results
write.table(res2, "correlation_annotated_and_expressed_miRNAs.tsv", col.names=T, row.names=T, sep="\t", quote=F)
write.table(res3, "correlation_known_interactions_only.tsv", col.names=T, row.names=T, sep="\t", quote=F)
