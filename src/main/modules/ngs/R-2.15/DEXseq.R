# TOOL DEXSeq.R: "Differential exon expression using DEXSeq" (Infers differential exon usage from RNA-seq data using the Bioconductor package DEXSeq. Replicates are necessary for this tool to work properly. In order to prepare the input, run the tool \"Count aligned reads per exons for DEXSeq\" for your BAM files and combine the results to a count table using the tool \"Utilities - Define NGS experiment\". Please use the group column of the phenodata file to indicate your experimental groups.)
# INPUT countfile.tsv: "Count table" TYPE GENERIC 
# INPUT phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT OPTIONAL DEXSeq-exons.pdf: DEXSeq-exons.pdf
# OUTPUT DEXSeq-result-table.tsv: DEXSeq-result-table.tsv
# PARAMETER organism: "Organism" TYPE [Homo_sapiens.GRCh37.68.chr.DEXSeq.gtf: "Human (hg19.68)", Mus_musculus.GRCm38.68.chr.DEXSeq.gtf: "Mouse (mm10.68)", Rattus_norvegicus.RGSC3.4.68.chr.DEXSeq.gtf: "Rat (rn4.68)"] DEFAULT Homo_sapiens.GRCh37.68.chr.DEXSeq.gtf (Which organism is your data from.)
# PARAMETER OPTIONAL dispersion: "Common dispersion" TYPE DECIMAL FROM 0 TO 100 DEFAULT 0.1 (This common dispersion value is used automatically for all exons, if dispersions can not be estimated. In this case no graphical output is generated.)

# JTT 15.7.2013


# Loads the library 
library(DEXSeq)

# Reads the phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
phenodata$condition<-phenodata$group
pd<-phenodata[,"group", drop=FALSE]
rownames(pd)<-phenodata$sample
colnames(pd)<-"condition"
pd$condition<-factor(pd$condition)

#if(any(as.vector(table(phenodata$group)))<2) {
#	stop("You need to have replicates for all groups you have specified.")
#}

# Path to the gff file
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
formuladispersion <- count ~ sample + condition * exon
ecs<-estimateDispersions(ecs, formula = formuladispersion)
ecs.fdf<-try(fitDispersionFunction(ecs), silent=TRUE)

if(class(ecs.fdf)=="try-error") {
   fData(ecs)$dispersion <- dispersion 
   doplot<-FALSE
} else {
   ecs<-ecs.fdf
   doplot<-TRUE
}

# Testing for differential exon usage
formula0<-count ~ sample + condition + exon
formula1<-count ~ sample + condition * I(exon==exonID) 
ecs <- testForDEU(ecs, formula0 = formula0, formula1 = formula1)
ecs <- estimatelog2FoldChanges(ecs)
res <- DEUresultTable(ecs)
write.table(res, "DEXSeq-result-table.tsv", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

# Visualization
if(doplot) {
   genes<-unique(as.character(res[which(res$padjust<=0.05),]$geneID))
   pdf("DEXSeq-exons.pdf",)
   for(i in 1:length(genes)) {
      plottry<-try(plotDEXSeq(ecs, genes[i], displayTranscripts = TRUE, cex.axis = 1.2, cex = 1.3, lwd = 2, legend = TRUE))
      if(class(plottry)=="try-error") {
         plot(x=1, y=1, xlab="", ylab="", axes=F, type="")
         title(main=genes[i])
         text(x=1, y=1, "No results to plot for this transcript")
      }
   }
   dev.off()
}
