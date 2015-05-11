# TOOL deseq2-normalized-counts-for-gene.R: "Plot normalized counts for a gene" (Plot normalized counts for a gene using the DESeq2 Bioconductor package.)
# INPUT data.tsv: "Count table" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata file" TYPE GENERIC
# OUTPUT OPTIONAL normalized_counts.pdf
# PARAMETER gene.name: "Gene name" TYPE STRING (Gene name.)
# PARAMETER OPTIONAL show.names: "Show names in plot" TYPE [yes, no] DEFAULT no (Show sample names in plot. In more cpmplex cases this may make the plot too cluttered.)

# AMS 21.4.2015 

# Loads the libraries
library(DESeq2)
library(ggplot2)


# Load the count table and extract expression value columns
dat <- read.table("data.tsv", header=T, sep="\t", row.names=1)
dat2 <- dat[,grep("chip", names(dat))]

# Get the experimental group information from the phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
condition <- as.character (phenodata[,pmatch("group",colnames(phenodata))])

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=dat2, colData=data.frame(condition), design = ~ condition)

dds <- DESeq(dds)
#res <- results(dds)
desc <- phenodata[,7]

# Make plot as pdf
pdf(file="normalized_counts.pdf")
	d <- plotCounts(dds, gene=gene.name, intgroup="condition",returnData=TRUE)
	if (show.names == "yes"){
		ggplot(d, aes(x=condition, y=count, log="y")) +
			#geom_point(color="blue", size=3, shape=5, position=position_jitter(w=0.1,h=0)) +
			geom_point(color="blue", size=3, shape=5) +
			geom_text(aes(label=desc),hjust=-0.5, vjust=0, color="black", size=4) +
			ylab("normalized counts") +
			xlab("group") +
			ggtitle(gene.name)	
	}else{
		ggplot(d, aes(x=condition, y=count, log="y")) +
			#geom_point(color="blue", size=3, shape=5, position=position_jitter(w=0.1,h=0)) +
			geom_point(color="blue", size=3, shape=5) +
			ylab("normalized counts") +
			xlab("group") +
			ggtitle(gene.name)
	}
dev.off()



