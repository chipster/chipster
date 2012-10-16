# TOOL DEXSeq.R: "Differential exon expression using DEXSeq" (Infers differential exon usage from RNA-seq data using the Bioconductor package DEXSeq. You can create the input count table and phenodata file using the tool Utilities - Define NGS experiment. In it highly recommended that all the groups defined by the phenodata group column have replicates.  If they do not, you have to estimate the dispersion manually by defining the parameter dispersion.)
# INPUT countfile.tsv: countfile.tsv TYPE GENERIC 
# INPUT phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT OPTIONAL DEXSeq-exons.pdf: DEXSeq-exons.pdf
# OUTPUT DEXSeq-result-table.tsv: DEXSeq-result-table.tsv
# PARAMETER organism: "Organism" TYPE [Homo_sapiens.GRCh37.68.chr.DEXSeq.gtf: "Human (hg19.68)", Mus_musculus.GRCm38.68.chr.DEXSeq.gtf: "Mouse (mm10.68)", Rattus_norvegicus.RGSC3.4.68.chr.DEXSeq.gtf: "Rat (rn4.68)"] DEFAULT Homo_sapiens.GRCh37.68.chr.DEXSeq.gtf (Which organism is your data from.)
# PARAMETER OPTIONAL dispersion: "Common dispersion" TYPE DECIMAL FROM 0 TO 100 DEFAULT 1 (The common dispersion used automatically for all transcripts, if there are no replicates, or if the correct way to estimate the dispersions runs into problem. If the dispersions can not be estimated, no graphical output is generated.)

# JTT 23.9.2012

# setwd("C:/Users/Jarno Tuimala/Desktop/dexseq")


# Loads the library 
library(DEXSeq)

# Reads the phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
phenodata$condition<-phenodata$group

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

ecs = read.HTSeqCounts(countfiles = dir(pattern="jtt"), design = phenodata, flattenedfile = gtf)
sampleNames(ecs)<-phenodata$original_name

# Normalization
ecs<-estimateSizeFactors(ecs)

# Estimate dispersion
formuladispersion <- count ~ sample + group * exon
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
formula0<-count ~ sample + group + exon
formula1<-count ~ sample + group * I(exon==exonID) 
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

