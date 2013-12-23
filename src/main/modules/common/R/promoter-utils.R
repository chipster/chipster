# Utility function for promoter retrieval. Used by the different promoter tools

retreive_promoters <- function(species, promoter.size, multiple.promoters, norm.file, pheno.file) {
	# Load Sequence library
	library(GenomicFeatures)

	# Loads the normalized data
	dat<-read.table(norm.file, header=T, sep="\t", row.names=1)
	phenodata<-read.table(pheno.file, header=T, sep="\t")

	# Loads the correct annotation library
	lib<-as.character(phenodata$chiptype[1])
	if (length(grep(".db", lib)) == 0 & length(grep("pmcdf", lib)) == 0) {
	   	lib <- paste(lib, ".db", sep="")
	}
	library(package=lib, character.only=T)

	# Creates a variable for environement
	if(species=="drosophila") {
	   	env<-paste(lib, "ENSEMBL", sep="")
	} else {
	   	env<-paste(lib, "ENTREZID", sep="")
	}
	env <- sub( ".db", "", env) # if chip contained ".db", remove it

	# Retrieve entrezid / ensembl list 
	refseq<-as.vector(unlist(mget(rownames(dat), envir=get(env))))
	#select(get(lib), rownames(dat), "ENTREZID", "PROBEID")

	if(species=="human") {
	   	genome <-  "BSgenome.Hsapiens.UCSC.hg19"
	   	database <- "TxDb.Hsapiens.UCSC.hg19.knownGene" 
	   	organism <- "Hsapiens"
	   	library(org.Hs.eg.db)
	}
	if(species=="mouse") {
	   	genome <-  "BSgenome.Mmusculus.UCSC.mm10"
	   	database <- "TxDb.Mmusculus.UCSC.mm10.knownGene" 
	   	organism <- "Mmusculus"
	   	library(org.Mm.eg.db)
	}
	if(species=="rat") {
	   	genome <-  "BSgenome.Rnorvegicus.UCSC.rn5"
	   	database <- "TxDb.Rnorvegicus.UCSC.rn5.refGene" 
	   	organism <- "Rnorvegicus"
	   	library(org.Rn.eg.db)
	}
	if(species=="drosophila") {
	   	genome <-  "BSgenome.Dmelanogaster.UCSC.dm3"
	   	database <- "TxDb.Dmelanogaster.UCSC.dm3.ensGene" 
	   	organism <- "Dmelanogaster"
	   	library(org.Dm.eg.db)
	}
	if(species=="yeast") {
	   	genome <-  "BSgenome.Scerevisiae.UCSC.sacCer3"
	   	database <- "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene" 
	   	organism <- "Scerevisiae"
		library(org.Sc.sgd.db)
	}

	# Extract Upstream
	library(database, character.only=T)
	library(genome, character.only=T)

	database.entries <- names(transcriptsBy(get(database), by="gene"))
	found.refseq <- refseq[!is.na(match(refseq, database.entries))]

	if(length(found.refseq) > 0) {
   		chromosomal.loc <- transcriptsBy(get(database), by="gene") [found.refseq]
   		promoter <- getPromoterSeq(chromosomal.loc, get(organism), upstream=as.integer(promoter.size), downstream=0)
	} else {
   		stop("CHIPSTER-NOTE: None of your input probes IDs were linked with Refseq IDs")
	}

	# Retrieving the sequences
	ret.values <- NULL;
	ret.values$seq <- NULL
	ret.values$seq.names <- NULL
	ret.values$nos.names <- NULL
	if(multiple.promoters=="no") {
	   	for(i in 1:length(promoter)) {
	      	ret.values$seq <- c(ret.values$seq, toString(promoter[[i]][1]))
	      	ret.values$seq.names <- c(ret.values$seq.names, paste(">", rownames(dat)[which(refseq==names(promoter)[i])] , " (Refseq ID: ", names(promoter)[i], ")", sep=""))
	   	}
	} else {
	   	for(i in 1:length(promoter)) {
	      	for(j in 1:length(promoter[[i]])) {
		      	ret.values$seq <- c(ret.values$seq, toString(promoter[[i]][j]))
		      	ret.values$seq.names <- c(ret.values$seq.names, paste(">", rownames(dat)[which(refseq==names(promoter)[i])] , ".", j, " (Refseq ID: ", names(promoter)[i], ")", sep=""))
	      	}
	   	}   
	}

	#print sequence information for probes without sequence information
	for(i in 1:length(rownames(dat))) {
	   	if(!is.na(match(refseq, database.entries))[i]==FALSE) {
	      	ret.values$nos.names <- c(ret.values$nos.names, paste(">", rownames(dat)[i], sep=""))
	   	}
	}

	return(ret.values)
}

