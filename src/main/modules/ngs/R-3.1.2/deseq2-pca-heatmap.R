# TOOL deseq2-pca-heatmap.R: "Experiment level QC for RNA-seq data" (Given a table of read counts for an experiment, this tool creates PCA and heatmap plots for the samples using the DESeq2 Bioconductor package. You can create the input count table and phenodata file using the tool \"Utilities - Define NGS experiment\".)
# INPUT data.tsv TYPE GENERIC
# INPUT phenodata.tsv TYPE GENERIC
# OUTPUT OPTIONAL pca-deseq2.bed
# OUTPUT OPTIONAL heatmap-deseq2.pdf
# PARAMETER column: "Column describing groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the experimental groups.)

# EK 3.2.2015 

# Loads the libraries
library(DESeq2)
library(gplots)
library(RColorBrewer)

# Load the count table and extract expression value columns
dat <- read.table("data.tsv", header=T, sep="\t", row.names=1)
dat2 <- dat[,grep("chip", names(dat))]

# Get the experimental group information from the phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
condition <- as.character (phenodata[,pmatch(column,colnames(phenodata))])

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=dat2, colData=data.frame(condition), design = ~ condition)

# Calculate size factors and dispersion values
dds <- estimateSizeFactors(dds) 
dds <- estimateDispersions(dds)	

# Perform transformation
vst<-varianceStabilizingTransformation(dds)
vstmat<-assay(vst)

# Make PCA plot as pdf
pdf(file="pca-deseq2.pdf")
plotPCA(vst,intgroup="condition")
dev.off()

# Make a distance matrix and name samples accroding to the phenodata description column
distvst<-dist(t(vstmat))
mdistvst<-as.matrix(distvst)
rownames(mdistvst)<-colnames(mdistvst)<-as.vector(phenodata$description)

# Calculate hierarchical clustering
hcvst<-hclust(distvst)

# Make the colors and plot the heatmap as pdf
hmcol<-colorRampPalette(brewer.pal(9,"GnBu"))(100)
pdf(file="heatmap-deseq2.pdf")
heatmap.2(mdistvst,Rowv=as.dendrogram(hcvst),Colv=as.dendrogram(hcvst),symm=TRUE,trace="none",col=rev(hmcol),margin=c(13,13))
dev.off()

# rld<-rlog(dds)
# rldmat<-assay(rld)
# distrld<-dist(t(rldmat))
# mdistrld<-as.matrix(distrld)
# hcrld<-hclust(distrld)




