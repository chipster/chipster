# ANALYSIS Utilities/"Search by gene name" (Search genes using gene symbol, rownames or chromosome locations. 
# Rowname must the name of the gene that appears
# on the rows of the data file. Chrosomome location must be a chromosome name, such as X. The "mode" option
# determines whether the queried gene name or chromosome should be included or excluded in the output data table.
# Note that the match when searching for a gene symbol can be exact, to get only that specific gene, or be a part of the
# match, to allow extraction of for example all members of a gene family (sharing the same part in their name.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT search.tsv
# PARAMETER search.for [rowname, symbol, chromosome.location] DEFAULT rowname (What to search with.)
# PARAMETER query STRING DEFAULT empty (Query word.)
# PARAMETER exact [yes, no] DEFAULT yes (Defines whether the search string should match exactly or be allowed to be a part of a match. Note, only applies for searches on symbol.)
# PARAMETER mode [include, exclude] DEFAULT include (Defines whether the found genes should be included or excluded from the resulting data table.)

# Search genes by rowname, symbol or chromosome location
# MG, 26.1.2011, based on a tool developed by JT, 4.7.2006

# Loads data (which file to search)
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1, quote="", comment.char="")

# Do the matching and write output
# For rownames
if (search.for == "rowname") {
	if(mode=="include") dat2<-dat[grep(query, row.names(dat), invert=FALSE),]
	if(mode=="exclude") dat2<-dat[grep(query, row.names(dat), invert=TRUE),]
	write.table(dat2, file=("search.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
}
# For symbols
if (search.for == "symbol") {
	# Make sure that the symbol column is in fact in the data
	if (length(grep("symbol", colnames(dat)))<1) {
		stop("CHIPSTER-NOTE: Your data does not have a column with information for symbols. For supported chip types please first run the Annotation / Add annotations to data tool to generate the necessary column.")
	}
	# Setup fuzzy or exact query
	query <- toupper(query)
	if (exact == "yes") {
		query <- paste("^", query, sep="")
	}
	if(mode=="include") dat2<-dat[grep(query, toupper(dat$symbol), invert=FALSE),]
	if(mode=="exclude") dat2<-dat[grep(query, toupper(dat$symbol), invert=TRUE),]
	write.table(dat2, file=("search.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
}
# For chromosome location
if(search.for == "chromosome.location") {
	# Loads libraries where applicable
	phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
	if(phenodata$chiptype[1]!="cDNA" & phenodata$chiptype[1]!="Illumina" & phenodata$chiptype[1]!="miRNA" & phenodata$chiptype[1]!="other" & phenodata$chiptype[1]!="empty") {
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
	# NOT supported for data of type "other"
	if (phenodata$chiptype[1]=="other" | phenodata$chiptype[1]=="empty") {
		stop("CHIPSTER-NOTE: To search by chromosome location is not supported for other than microarray data. Please rerun the tool with the search.for parameter set to rownames or symbol.")
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

# EOF