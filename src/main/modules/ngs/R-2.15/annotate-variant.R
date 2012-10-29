# TOOL annotate-variant.R: "Annotate variants" (Annotate variants listed in a VCF file using the Bioconductor package VariantAnotation.)
# INPUT input.vcf: "Sorted or unsorted VCF file" TYPE GENERIC
# OUTPUT all-variants.tsv
# OUTPUT coding-variants.tsv
# PARAMETER genome: "Genome" TYPE [hg19: "Human (hg19\)"] DEFAULT hg19 (Reference sequence)


# 31.8.2012 JTT

#input.vcf<-system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#input.vcf<-"var.flt2.vcf"
#genome<-"hg19"
#setwd("C:/Users/Jarno Tuimala/Desktop/c")


# Read data
library(VariantAnnotation)
vcf<-readVcf("input.vcf", genome)

# Correct the chromosome names 
#vcf@rowData@seqnames@values
if(length(grep("chr", vcf@rowData@seqnames@values))>=1) {
   vcf2<-vcf
} else {
   #vcf@rowData@seqnames@values<-paste("chr", vcf@rowData@seqnames@values, sep="")
   vcf2<-vcf
   if(genome=="hg19") {
      vcf2 <- renameSeqlevels(vcf, c("1"="chr1", "2"="chr2", "3"="chr3", "4"="chr4", "5"="chr5", "6"="chr6", "7"="chr7", "8"="chr8", "9"="chr9", "10"="chr10",
                                     "11"="chr11", "12"="chr12", "13"="chr13", "14"="chr14", "15"="chr15", "16"="chr16", "17"="chr17", "18"="chr18", "19"="chr19", "20"="chr20",
                                     "21"="chr21", "22"="chr22", "Y"="chrY", "X"="chrX"))
   }
}

# Exon isoform database
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Locate coding variant
# loc <- locateVariants(rd, txdb, CodingVariants())

# Locate splicing variants
# spl<-locateVariants(rd, txdb, SpliceSiteVariants())

# Locate all variants
rd<-rowData(vcf2)
allvar <- locateVariants(rd, txdb, AllVariants())
allvar2<-as.data.frame(allvar)
colnames(allvar2)<-toupper(colnames(allvar2))

# Get gene annotations from EntrezIDs
ig<-allvar2[!is.na(allvar2$GENEID),]
nig<-allvar2[is.na(allvar2$GENEID),]
library(org.Hs.eg.db)
symbol <- select(org.Hs.eg.db, keys=unique(ig[ig$LOCATION=="coding",]$GENEID), keytype="ENTREZID", cols="SYMBOL")
genename <- select(org.Hs.eg.db, keys=unique(ig[ig$LOCATION=="coding",]$GENEID), keytype="ENTREZID", cols="GENENAME")
ensg <- select(org.Hs.eg.db, keys=unique(ig[ig$LOCATION=="coding",]$GENEID), keytype="ENTREZID", cols="ENSEMBL")
ig2<-merge(ig, symbol, by.x="GENEID", by.y="ENTREZID", all.x=TRUE)
ig3<-merge(ig2, genename, by.x="GENEID", by.y="ENTREZID", all.x=TRUE)
ig4<-merge(ig3, ensg, by.x="GENEID", by.y="ENTREZID", all.x=TRUE)
nig$SYMBOL<-rep(NA, nrow(nig))
nig$GENENAME<-rep(NA, nrow(nig))
nig$ENSEMBL<-rep(NA, nrow(nig))
allvar3<-rbind(ig4, nig)

# Predict coding amino acid changes
library(BSgenome.Hsapiens.UCSC.hg19)
coding <- predictCoding(vcf2, txdb, seqSource=Hsapiens)
cod<-elementMetadata(coding)
cod2<-as.list(cod)
names(cod2)<-toupper(names(cod2))
#cod3<-data.frame(cod2$GENEID, cod2$CDSID, cod2$TXID, cod2$CONSEQUENCE, as.data.frame(cod2$PROTEINLOC), as.data.frame(cod2$CDSLOC)[1], as.data.frame(cod2$CDSLOC)[2], as.data.frame(cod2$CDSLOC)[3], as.data.frame(cod2$varAllele), as.data.frame(cod2$REFCODON), as.data.frame(cod2$VARCODON))
cod3<-data.frame(geneID=cod2$GENEID, cdsID=cod2$CDSID, txID=cod2$TXID, consequence=cod2$CONSEQUENCE, cdsStart=as.data.frame(cod2$CDSLOC)[,1], cdsEnd=as.data.frame(cod2$CDSLOC)[,2], width=as.data.frame(cod2$CDSLOC)[,3], varAllele=as.data.frame(cod2$VARALLELE)[,1], refCodon=as.data.frame(cod2$REFCODON)[,1], varCodon=as.data.frame(cod2$VARCODON)[,1], refAA=as.data.frame(cod2$REFAA)[,1], varAA=as.data.frame(cod2$VARAA)[,1])
symbol <- select(org.Hs.eg.db, keys=unique(cod3$geneID), keytype="ENTREZID", cols="SYMBOL")
genename <- select(org.Hs.eg.db, keys=unique(cod3$geneID), keytype="ENTREZID", cols="GENENAME")
ensg <- select(org.Hs.eg.db, keys=unique(cod3$geneID), keytype="ENTREZID", cols="ENSEMBL")
cod32<-merge(cod3, symbol, by.x="geneID", by.y="ENTREZID", all.x=TRUE)
cod33<-merge(cod32, genename, by.x="geneID", by.y="ENTREZID", all.x=TRUE)
cod34<-merge(cod33, ensg, by.x="geneID", by.y="ENTREZID", all.x=TRUE)

# Write results to disk
write.table(allvar3, "all-variants.tsv", col.names=T, row.names=F, sep="\t", quote=FALSE)
write.table(cod34, "coding-variants.tsv", col.names=T, row.names=F, sep="\t", quote=FALSE)
