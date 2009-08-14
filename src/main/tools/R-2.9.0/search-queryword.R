# ANALYSIS Utilities/"Search by gene name" (Search genes using gene names, rownames or chromosome locations. 
# Gene name must the gene's HUGO name, such as XRCC1. Rowname must the name of the gene that appears
# on the rows of all data files. Chrosomome location must be a chromosome name, such as X.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT search.tsv
# PARAMETER search.for [rowname, Genename, ChromosomeLocation] DEFAULT rowname (What to search with)
# PARAMETER query STRING DEFAULT empty (Query word)


# Search genes by name, AffyID, correlation or chromosome location
# JTT 4.7.2006

# Renaming variables
meth<-search.for
query<-query

# Loads libraries
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
if(phenodata$chiptype[1]!="cDNA" | phenodata$chiptype[1]!="Illumina") {
   lib<-phenodata$chiptype[1]
   lib<-as.character(lib)
}

# Loads the selected chip file
if(phenodata$chiptype[1]!="cDNA" | phenodata$chiptype[1]!="Illumina") {
   library(package=lib, character.only=T)
}

# Loads data (which file to search)
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

if(meth=="rowname") {
   dat2<-dat[grep(query, row.names(dat)),]
   write.table(dat2, file=("search.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
}

if(meth=="Genename") {
   query<-paste("^", query, sep="")
   lib2<-sub('.db','',lib)
   env<-paste(lib2, "SYMBOL", sep="")
   env2<-get(env)
   list<-unlist(as.list(env2))[grep(query, unlist(as.list(env2)))]
   affyid<-names(list)
   len<-length(affyid)
   dat2<-c()
   for(i in 1:len) {
      dat2<-rbind(dat2, dat[which(row.names(dat)==affyid[i]),])
   }
   write.table(dat2, file=("search.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
}

if(meth=="ChromosomeLocation") {
   lib2<-sub('.db','',lib)
   env<-paste(lib2, "CHR", sep="")
   env2<-get(env)
   affyids<-data.frame(id=names(as.list(env2))[grep(query, as.list(env2))])
   dat2<-merge(dat, affyids, by.x="row.names", by.y="id")
   row.names(dat2)<-dat2$Row.names
   dat2<-dat2[,-1]
   write.table(dat2, file="search.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}




