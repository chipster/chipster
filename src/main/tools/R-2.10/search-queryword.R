# ANALYSIS Utilities/"Search by gene name" (Search genes using gene names, rownames or chromosome locations. 
# Gene name must be the gene's HUGO name, such as XRCC1. Rowname must the name of the gene that appears
# on the rows of all data files. Chrosomome location must be a chromosome name, such as X. The "mode" option
# determines whether the queried gene name or chromosome should be included or excluded in the output data table.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT search.tsv
# PARAMETER search.for [rowname, Genename, ChromosomeLocation] DEFAULT rowname (What to search with)
# PARAMETER query STRING DEFAULT empty (Query word)
# PARAMETER mode [include, exclude] DEFAULT include (Defines whether the found genes should be included or excluded from the resulting data table.)

# Search genes by name, AffyID, correlation or chromosome location
# JTT 4.7.2006
#
# MG, 23.2.2010, modified to allow option to exclude query genes or chromosomes
#
# MG, 27.12.2010 modifed to cope with other than microarray data
#

# Renaming variables
meth <- search.for
query <- query

# Loads libraries where applicable
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
if(phenodata$chiptype[1]!="cDNA" & phenodata$chiptype[1]!="Illumina" & phenodata$chiptype[1]!="miRNA" & phenodata$chiptype[1]!="other" & phenodata$chiptype=="empty") {
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
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

if(meth=="rowname") {
	if(mode=="include") dat2<-dat[grep(query, row.names(dat), invert=FALSE),]
	if(mode=="exclude") dat2<-dat[grep(query, row.names(dat), invert=TRUE),]
	write.table(dat2, file=("search.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
}

if(meth=="Genename") {
	# NOT supported for data of type "other"
	if (phenodata$chiptype[1]=="other" | phenodata$chiptype[1]=="empty") {
		stop("CHIPSTER-NOTE: To search by gene name is not supported for other than microarray data. Please rerun the tool with the rownames parameter enabled.")
	}
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
	# NOT supported for data of type "other"
	if (phenodata$chiptype[1]=="other" | phenodata$chiptype[1]=="empty") {
		stop("CHIPSTER-NOTE: To search by gene name is not supported for other than microarray data. Please rerun the tool with the rownames parameter enabled.")
	}
   lib2<-sub('.db','',lib)
   env<-paste(lib2, "CHR", sep="")
   env2<-get(env)
   if(mode=="include") affyids<-data.frame(id=names(as.list(env2))[grep(query, as.list(env2), invert=FALSE)])
   if(mode=="exclude") affyids<-data.frame(id=names(as.list(env2))[grep(query, as.list(env2), invert=TRUE)])
      dat2<-merge(dat, affyids, by.x="row.names", by.y="id")
   row.names(dat2)<-dat2$Row.names
   dat2<-dat2[,-1]
   write.table(dat2, file="search.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}




