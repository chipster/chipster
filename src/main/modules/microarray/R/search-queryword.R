# TOOL search-queryword.R: "Search by gene name" (Search genes using gene names, rownames or chromosome locations. Gene name must be the gene's HUGO name, such as XRCC1. Rowname must the name of the gene that appears on the rows of all data files. Chrosomome location must be a chromosome name, such as X or 1. The mode option determines whether the queried gene name or chromosome should be included or excluded in the output data table.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT search.tsv: search.tsv 
# PARAMETER search.for: "Search for" TYPE [rowname: rowname, symbol: symbol, description: description, Genename: "Gene name", ChromosomeLocation: "Chromosome Location"] DEFAULT rowname (What to search with)
# PARAMETER query: Query TYPE STRING DEFAULT empty (Query word)
# PARAMETER mode: Mode TYPE [include: include, exclude: exclude] DEFAULT include (Defines whether the found genes should be included or excluded from the resulting data table.)

# JTT, 4.7.2006: Search genes by name, AffyID, correlation or chromosome location
# MG, 23.2.2010: to allow option to exclude query genes or chromosomes
# IS, 12.10.2012: to cope with tables with gene descriptions (that typically contain 's)

# Renaming variables
meth<-search.for
query<-query

# Loads libraries
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
if(phenodata$chiptype[1]!="cDNA" & phenodata$chiptype[1]!="Illumina" & phenodata$chiptype[1]!="miRNA" ) {
   lib<-phenodata$chiptype[1]
   lib<-as.character(lib)
   
   # Account for the fact that annotation packages are from version 2.3 of Bioconductor
   # named with an ".db" suffix. Add the suffix when missing to support data files
   # from Chipster 1.3 and earlier. 
   if (length(grep(".db", lib)) == 0 & length(grep("pmcdf", lib)) == 0) {
        lib <- paste(lib, ".db", sep="")
   }
   
   library(package=lib, character.only=T)
}

# Loads data (which file to search)
file <- 'normalized.tsv'
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, check.names=FALSE)

options(scipen=10)

if(meth=="rowname") {
	if(mode=="include") dat2<-dat[grep(query, row.names(dat), invert=FALSE),]
	if(mode=="exclude") dat2<-dat[grep(query, row.names(dat), invert=TRUE),]
	write.table(dat2, file=("search.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
}

if(meth=="symbol") {
   if(length(grep("symbol", colnames(dat)))==1) {
      if(mode=="include") dat2<-dat[grep(query, dat[,"symbol"], invert=FALSE),]
      if(mode=="exclude") dat2<-dat[grep(query, dat[,"symbol"], invert=TRUE),]
      write.table(dat2, file=("search.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
   } else {
      stop("CHIPSTER-NOTE: Your data does not have symbol filed")
   }
}

if(meth=="description") {
   if(length(grep("description", colnames(dat)))==1) {
      if(mode=="include") dat2<-dat[grep(query, dat[,"description"], invert=FALSE),]
      if(mode=="exclude") dat2<-dat[grep(query, dat[,"description"], invert=TRUE),]
      write.table(dat2, file=("search.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
   } else {
      stop("CHIPSTER-NOTE: Your data does not have description filed")
   }
}


if(meth=="Genename") {
   query<-paste("^", query, sep="")
   lib2<-sub('.db','',lib)
   env<-paste(lib2, "SYMBOL", sep="")
   env2<-get(env)
   list<-unlist(as.list(env2))[grep(query, unlist(as.list(env2)), invert=FALSE)]
	affyid<-names(list)
	matching_indices <- match(affyid, rownames(dat))
	if (mode=="include") dat2 <- dat[matching_indices,]
	if (mode=="exclude") dat2 <- dat[-matching_indices,]
   write.table(dat2, file=("search.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
}

if(meth=="ChromosomeLocation") {
   lib2<-sub('.db','',lib)
   env<-paste(lib2, "CHR", sep="")
   #env<-paste(lib2, "CHRLOC", sep="")
   env2<-get(env)
   if(mode=="include") affyids<-data.frame(id=names(as.list(env2))[grep(query, as.list(env2), invert=FALSE)])
   if(mode=="exclude") affyids<-data.frame(id=names(as.list(env2))[grep(query, as.list(env2), invert=TRUE)])
      dat2<-merge(dat, affyids, by.x="row.names", by.y="id")
   row.names(dat2)<-dat2$Row.names
   dat2<-dat2[,-1]
   write.table(dat2, file="search.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# EOF
