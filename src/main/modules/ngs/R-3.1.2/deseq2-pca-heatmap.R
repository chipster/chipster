# TOOL deseq2-pca-heatmap.R: "PCA and heatmap of samples with DESeq2" (Creates PCA and heatmap plots for samples using the DESeq2 Bioconductor package. Visualizing similarities and dissimilarities between samples can help you to perform experiment level quality control. This tool takes as input a table of raw counts. The count table has to be associated with a phenodata file describing the experimental groups. You can create the count table and phenodata file using the tool \"Utilities - Define NGS experiment\".)
# INPUT data.tsv: "Count table" TYPE GENERIC
# INPUT phenodata.tsv: "Phenodata file" TYPE GENERIC
# OUTPUT OPTIONAL pca-deseq2.pdf
# OUTPUT OPTIONAL heatmap-deseq2.pdf
# PARAMETER column: "Column describing groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the experimental groups.)

# EK 3.2.2015 

# Loads the libraries
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(ggplot2)

# Load the count table and extract expression value columns
dat <- read.table("data.tsv", header=T, sep="\t", row.names=1)
dat2 <- dat[,grep("chip", names(dat))]

# Get the experimental group information from the phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
condition <- as.character (phenodata[,pmatch(column,colnames(phenodata))])

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=dat2, colData=data.frame(condition), design = ~ condition)

# Calculate size factors and estimate dispersions
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# Perform transformation
vst<-varianceStabilizingTransformation(dds)
vstmat<-assay(vst)

# Make PCA plot as pdf
# plotPCA(vst,intgroup="condition")

data <- plotPCA(vst, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pdf(file="pca-deseq2.pdf")
ggplot(data, aes(PC1, PC2, color=condition)) +
		geom_point(size=6,shape=0) +
		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
		ylab(paste0("PC2: ",percentVar[2],"% variance"))
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

# To be added: rlog transformation
# rld<-rlog(dds)
# rldmat<-assay(rld)
# distrld<-dist(t(rldmat))
# mdistrld<-as.matrix(distrld)
# hcrld<-hclust(distrld)




