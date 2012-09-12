# TOOL DEXSeq.R: "Differential exon usage in RNA-seq" (Infer differential exon usage from RNA-seq data. Input is an NGS experiment defined by the tool in the utilities category.)
# INPUT countfile.tsv: countfile.tsv TYPE GENERIC 
# INPUT phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT DEXSeq-exons.pdf: DEXSeq-exons.pdf
# OUTPUT DEXSeq-result-table.tsv: DEXSeq-result-table.tsv
# PARAMETER organism: "Organism" TYPE [Homo_sapiens.GRCh37.62.chr.gtf: "Human (hg19.62)", Mus_musculus.NCBIM37.62.chr.gtf: "Mouse (mm9.62)", Rattus_norvegicus.RGSC3.4.62.chr.gtf: "Rat (rn4.62)"] DEFAULT Homo_sapiens.GRCh37.62.chr.gtf (Which organism is your data from.)


# JTT 14.8.2012

# Loads the library 
library(DEXSeq)

# Reads the phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")

# Path to the gff file
gtf <- file.path(chipster.tools.path, "genomes", organism)

# Reads the data
d<-read.table("countfile.tsv", header=TRUE, sep="\t")
d2<-d[,grep("chip", colnames(d))]
cn<-substr(colnames(d2), 6, nchar(colnames(d2)))
for(i in 1:ncol(d2)) {
   v<-d2[,i, drop=F]
   rownames(v)<-rownames(d2)
   write.table(v, paste(cn[i], ".jtt", sep=""), col.names=FALSE, row.names=TRUE, sep="\t", quote=FALSE)
}

ecs = read.HTSeqCounts(countfiles = dir(pattern="jtt"), design = phenodata, flattenedfile = gtf)
sampleNames(ecs)<-phenodata$original_name

# Normalization
ecs<-estimateSizeFactors(ecs)

# Estimate dispersion
phenodata$condition<-phenodata$group
formuladispersion <- count ~ sample + group * exon
ecs<-estimateDispersions(ecs, formula = formuladispersion)
ecs<-fitDispersionFunction(ecs)

# Testing for differential exon usage
formula0<-count ~ sample + group + exon
formula1<-count ~ sample + group * I(exon==exonID) 
ecs <- testForDEU(ecs, formula0 = formula0, formula1 = formula1)
ecs <- estimatelog2FoldChanges(ecs)
res <- DEUresultTable(ecs)
write.table(res, "DEXSeq-result-table.tsv", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

# Visualization
genes<-unique(as.character(na.omit(res[res$padjust<=0.05,])$geneID))

pdf("DEXSeq-exons.pdf",)
for(i in 1:length(genes)) {
   plotDEXSeq(ecs, genes[i], displayTranscripts = TRUE, cex.axis = 1.2, cex = 1.3, lwd = 2, legend = TRUE)
}
dev.off()
