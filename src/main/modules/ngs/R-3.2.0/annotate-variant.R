# TOOL annotate-variant.R: "Annotate variants" (Annotate variants listed in a VCF file using the Bioconductor package VariantAnnotation. Coding, UTR, intronic and promoter variants are annotated, intergenic variants are ignored.)
# INPUT input.vcf: "Sorted or unsorted VCF file" TYPE GENERIC
# OUTPUT all-variants.tsv
# OUTPUT coding-variants.tsv
# PARAMETER genome: "Genome" TYPE [hg19: "Human (hg19\)"] DEFAULT hg19 (Reference sequence)


# 31.8.2012 	JTT
# 15.11.2012 	JTT,EK  	Uses VariantAnnotation package 1.4.3. 
# 01.06.2015 	ML 			Modifications to move the tool to new R version
# 20.7.2015 	ML 	testaa 3


# Read data
library(VariantAnnotation)
vcf<-readVcf("input.vcf", genome)
vcf@rowData@seqnames@values <- factor(vcf@rowData@seqnames@values)  					#hmmm...... pidä kommentoituna.

# Correct the chromosome names: 
if(length(grep("chr", vcf@rowData@seqnames@values))>=1) {
   vcf2<-vcf
   rd<-rowRanges(vcf)
} else {
   if(genome=="hg19") {
		vcf2 <- vcf
		seqlevelsStyle(vcf2) <- "UCSC" 

      	rd<-rowRanges(vcf)
		rd2 <- rd
		seqlevelsStyle(rd2) <- "UCSC" 
   }
}

# Exon isoform database
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# remove mitochonrial DNA (it is causing problems with genome versions)
# vcf2 <- keepSeqlevels(vcf2, seqlevels(vcf2)[1:24])								# hmmmm.... ongelmia jos järjestys väärä.
vcf2 <- keepSeqlevels(vcf2, c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))



# Locate coding variant
# loc <- locateVariants(rd, txdb, CodingVariants())

# Locate splicing variants
# spl<-locateVariants(rd, txdb, SpliceSiteVariants())

# Locate all variants
codvar <- locateVariants(vcf2, txdb, CodingVariants())
intvar <- locateVariants(vcf2, txdb, IntronVariants())
utr5var <- locateVariants(vcf2, txdb, FiveUTRVariants())
utr3var <- locateVariants(vcf2, txdb, ThreeUTRVariants())
#intgvar <- locateVariants(vcf2, txdb, IntergenicVariants())
spsgvar <- locateVariants(vcf2, txdb, SpliceSiteVariants())
promvar <- locateVariants(vcf2, txdb, PromoterVariants())
allvar<-c(codvar, intvar, utr5var, utr3var, spsgvar, promvar)
#allvar <- locateVariants(rd, txdb, AllVariants())

allvar2<-as.data.frame(allvar, row.names=1:length(allvar))
colnames(allvar2)<-toupper(colnames(allvar2))

# Get gene annotations from EntrezIDs
ig<-allvar2[!is.na(allvar2$GENEID),]
nig<-allvar2[is.na(allvar2$GENEID),]
library(org.Hs.eg.db)
symbol <- select(org.Hs.eg.db, keys=unique(ig[ig$LOCATION=="coding",]$GENEID), keytype="ENTREZID", columns="SYMBOL")
genename <- select(org.Hs.eg.db, keys=unique(ig[ig$LOCATION=="coding",]$GENEID), keytype="ENTREZID", columns="GENENAME")
ensg <- select(org.Hs.eg.db, keys=unique(ig[ig$LOCATION=="coding",]$GENEID), keytype="ENTREZID", columns="ENSEMBL")
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
# cod3<-data.frame(geneID=cod2$GENEID, cdsID=cod2$CDSID, txID=cod2$TXID, consequence=cod2$CONSEQUENCE, cdsStart=as.data.frame(cod2$CDSLOC)[,1], cdsEnd=as.data.frame(cod2$CDSLOC)[,2], width=as.data.frame(cod2$CDSLOC)[,3], varAllele=as.data.frame(cod2$VARALLELE)[,1], refCodon=as.data.frame(cod2$REFCODON)[,1], varCodon=as.data.frame(cod2$VARCODON)[,1], refAA=as.data.frame(cod2$REFAA)[,1], varAA=as.data.frame(cod2$VARAA)[,1])
cod3<-data.frame(geneID=cod2$GENEID, cdsID=sapply(cod2$CDSID, FUN=function(x) paste(x, collapse=", ")), txID=cod2$TXID, consequence=cod2$CONSEQUENCE, cdsStart=as.data.frame(cod2$CDSLOC)[,1], cdsEnd=as.data.frame(cod2$CDSLOC)[,2], width=as.data.frame(cod2$CDSLOC)[,3], varAllele=as.data.frame(cod2$VARALLELE)[,1], refCodon=as.data.frame(cod2$REFCODON)[,1], varCodon=as.data.frame(cod2$VARCODON)[,1], refAA=as.data.frame(cod2$REFAA)[,1], varAA=as.data.frame(cod2$VARAA)[,1])


symbol <- select(org.Hs.eg.db, keys=as.character(unique(cod3$geneID)), keytype="ENTREZID", columns="SYMBOL")
genename <- select(org.Hs.eg.db, keys=as.character(unique(cod3$geneID)), keytype="ENTREZID", columns="GENENAME")
ensg <- select(org.Hs.eg.db, keys=as.character(unique(cod3$geneID)), keytype="ENTREZID", columns="ENSEMBL")
cod32<-merge(cod3, symbol, by.x="geneID", by.y="ENTREZID", all.x=TRUE)
cod33<-merge(cod32, genename, by.x="geneID", by.y="ENTREZID", all.x=TRUE)
cod34<-merge(cod33, ensg, by.x="geneID", by.y="ENTREZID", all.x=TRUE)

# Write results to disk
write.table(allvar3, "all-variants.tsv", col.names=T, row.names=F, sep="\t", quote=FALSE)
write.table(cod34, "coding-variants.tsv", col.names=T, row.names=F, sep="\t", quote=FALSE)