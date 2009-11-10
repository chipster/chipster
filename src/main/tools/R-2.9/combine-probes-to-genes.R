# ANALYSIS Utilities/"Combine probes to genes" (Calculates an average for probes or probesets for each gene in the
# dataset. The data file has to have a symbol column for this to work correctly. After running this tool, only
# expression values and gene symbols are retained in the data, all other columns and information are lost.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT combined.tsv
# PARAMETER produce.identifiers [no, yes] DEFAULT no (Should approximate identifiers for the gene be returned)


# Combine probes or probeset to genes
# 18.12.2008 JTT


# Loads the file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values from other data
dat2<-dat[,grep("chip", names(dat))]

# Preparations
test1<-aggregate(dat2[,1], list(dat$symbol), mean)
dat3<-matrix(ncol=ncol(dat2), nrow=nrow(test1), data=NA)
# dat4<-matrix(ncol=ncol(dat2), nrow=nrow(test1), data=NA)

# Combination
for(i in 1:ncol(dat2)) {
   m<-aggregate(dat2[,i], list(dat$symbol), mean)
   # s<-aggregate(dat2[,i], list(dat$symbol), sd)
   # s$x[is.na(s$x)]<-0
   dat3[,i]<-m$x
   # dat4[,i]<-s$x
}

# Second round of combination
if(produce.identifiers=="yes") {

   symbol<-m$Group.1

   # Generating rownames
   genes<-rep(NA, length(symbol)) 
   for(i in 1:length(symbol)) {
      genes[i]<-rownames(dat)[grep(symbol[i], dat$symbol)][1]
   }

   # Second round of combination, now according to rownames
   test2<-aggregate(dat3[,1], list(genes), mean)
   dat6<-matrix(ncol=ncol(dat2), nrow=nrow(test2), data=NA)

   for(i in 1:ncol(dat2)) {
      m<-aggregate(dat3[,i], list(genes), mean)
      dat6[,i]<-m$x
   }
   
   rownames(dat6)<-genes[!duplicated(genes)]
   colnames(dat6)<-colnames(dat2)

   write.table(data.frame(dat6), file="combined.tsv", col.names=T, quote=F, sep="\t", row.names=T)

} else {
 
   # Putting a data file together
   colnames(dat3)<-colnames(dat2)
   # colnames(dat4)<-paste("sd.", colnames(dat2), sep="")
   symbol<-m$Group.1
   # dat5<-data.frame(symbol, dat3, dat4) 
   dat5<-data.frame(dat3) 
   rownames(dat5)<-symbol

   # Writing data to disk
   write.table(data.frame(dat5), file="combined.tsv", col.names=T, quote=F, sep="\t", row.names=T)
}
