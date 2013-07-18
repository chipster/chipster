# TOOL DEXSeq.R: "Differential exon expression using DEXSeq" (Infers differential exon usage from RNA-seq data using the Bioconductor package DEXSeq. Replicates are necessary for this tool to work properly. In order to prepare the input, run the tool \"Count aligned reads per exons for DEXSeq\" for your BAM files and combine the results to a count table using the tool \"Utilities - Define NGS experiment\". Please use the group column of the phenodata file to indicate your experimental groups.)
# INPUT countfile.tsv: "Count table" TYPE GENERIC 
# INPUT phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT dexseq-all-genes.tsv: dexseq-all-genes.tsv
# OUTPUT OPTIONAL dexseq-genes-with-significant-exons.tsv: dexseq-genes-with-significant-exons.tsv
# OUTPUT OPTIONAL dexseq-exons.pdf: dexseq-exons.pdf
# OUTPUT OPTIONAL dexseq-MAplot.pdf: dexseq-MAplot.pdf
# OUTPUT OPTIONAL dexseq-dispersion-plot.pdf: dexseq-dispersion-plot.pdf
# PARAMETER organism: "Organism" TYPE [Homo_sapiens.GRCh37.68.chr.DEXSeq.gtf: "Human (hg19.68)", Mus_musculus.GRCm38.68.chr.DEXSeq.gtf: "Mouse (mm10.68)", Rattus_norvegicus.RGSC3.4.68.chr.DEXSeq.gtf: "Rat (rn4.68)"] DEFAULT Homo_sapiens.GRCh37.68.chr.DEXSeq.gtf (Which organism is your data from.)
# PARAMETER pvalue: "Threshold for adjusted p-value" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Threshold for BH adjusted p-values. If a gene has at least one exon below this p-value, all its exons will be included in the result list.)
# PARAMETER OPTIONAL dispersion: "Common dispersion" TYPE DECIMAL FROM 0 TO 100 DEFAULT 0.1 (If dispersions can not be estimated, this common dispersion value is used for all exons. In this case no graphical output is generated.)


# JTT 18.7.2013

####
####setwd("C:/Users/Jarno Tuimala/Desktop/DEXSeq")
####pvalue<-0.05
####dispersion<-0.1
####

# Loads the library 
library(DEXSeq)

# Reads the phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
phenodata$condition<-phenodata$group
pd<-phenodata[,"group", drop=FALSE]
rownames(pd)<-phenodata$sample
colnames(pd)<-"condition"
pd$condition<-factor(pd$condition)

if(any(as.vector(table(phenodata$group))<2)) {
	stop("You need to have replicates for all groups you have specified.")
}

# Path to the gff file
####
####gtf<-"C:/Users/Jarno Tuimala/Desktop/DEXSeq/Homo_sapiens.GRCh37.68.chr.DEXSeq.gtf"
####
gtf <- file.path(chipster.tools.path, "genomes", "gtf", organism)

# Reads the data
d<-read.table("countfile.tsv", header=TRUE, sep="\t")
d2<-d[,grep("chip", colnames(d))]
cn<-substr(colnames(d2), 6, nchar(colnames(d2)))
for(i in 1:ncol(d2)) {
   v<-d2[,i, drop=F]
   rownames(v)<-rownames(d2)
   write.table(v, paste(cn[i], ".jtt", sep=""), col.names=FALSE, row.names=TRUE, sep="\t", quote=FALSE)
}

ecs = read.HTSeqCounts(countfiles = dir(pattern="jtt"), design = pd, flattenedfile = gtf)
sampleNames(ecs)<-phenodata$original_name

# Normalization
ecs<-estimateSizeFactors(ecs)

# Estimate dispersion
####formuladispersion <- count ~ sample + condition * exon
####ecs<-estimateDispersions(ecs, formula = formuladispersion)
####ecs.fdf<-try(fitDispersionFunction(ecs), silent=TRUE)
ecs<-estimateDispersions(ecs)
ecs.fdf<-try(fitDispersionFunction(ecs), silent=TRUE)

if(class(ecs.fdf)=="try-error") {
   fData(ecs)$dispersion <- dispersion 
   doplot<-FALSE
} else {
   ecs<-ecs.fdf
   doplot<-TRUE
}

# Testing for differential exon usage
####formula0<-count ~ sample + condition + exon
####formula1<-count ~ sample + condition * I(exon==exonID) 
ecs <- testForDEU(ecs)
ecs <- estimatelog2FoldChanges(ecs)
res <- DEUresultTable(ecs)

siggenes<-as.character(unique(res$geneID[res$padjust<pvalue]))
res2<-res[as.character(res$geneID) %in% siggenes,]

write.table(res, "dexseq-all-genes.tsv", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
write.table(res2, "dexseq-genes-with-significant-exons.tsv", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

# Visualization
if(doplot & nrow(res2)>0) {
   genes<-unique(as.character(res[which(res$padjust<=pvalue),]$geneID))
   pdf("dexseq-exons.pdf", width=297/25.4, height=210/25.4)
   for(i in 1:length(genes)) {
      plottry<-try(plotDEXSeq(ecs, genes[i], displayTranscripts = FALSE, cex.axis = 1.2, cex = 1.3, lwd = 2, legend = TRUE))
      if(class(plottry)=="try-error") {
         plot(x=1, y=1, xlab="", ylab="", axes=F, type="")
         title(main=genes[i])
         text(x=1, y=1, "No results to plot for this gene")
      }
   }
   dev.off()
}

if(nrow(res)>0) {
   pdf("dexseq-MAplot.pdf", width=297/25.4, height=210/25.4)
   plot(x=log2(res[,6]), y=res[,7], xlab="mean of normalized counts (log2)", ylab="log2 fold change", type="n")
   points(x=log2(res[res[,5]>pvalue,6]), y=res[res[,5]>pvalue,7], col="black", pch=16, cex=0.5)
   points(x=log2(res[res[,5]<=pvalue,6]), y=res[res[,5]<=pvalue,7], col="#CC0000", pch=16, cex=0.5)
   dev.off()
}

plotDispEsts = function( cds, ymin, linecol="#ff000080",
  xlab = "mean of normalized counts", ylab = "dispersion",
  log = "xy", cex = 0.45, ... )
{
  px = rowMeans( counts( cds, normalized=TRUE ) )
  sel = (px>0)
  px = px[sel]

  py = fData(cds)$dispBeforeSharing[sel]
  if(missing(ymin))
      ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)

  plot(px, pmax(py, ymin), xlab=xlab, ylab=ylab,
    log=log, pch=ifelse(py<ymin, 6, 16), cex=cex, ... )
  xg = 10^seq( -.5, 5, length.out=100 )
  fun = function(x) { cds@dispFitCoefs[1] + cds@dispFitCoefs[2] / x }
  lines( xg, fun(xg), col=linecol, lwd=4)
}

pdf("dexseq-dispersion-plot.pdf", width=297/25.4, height=210/25.4)
plotDispEsts(ecs)
dev.off()
