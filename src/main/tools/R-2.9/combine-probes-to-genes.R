# ANALYSIS Utilities/"Combine probes to genes" (Calculates an average for probes or probesets for each gene in the
# dataset. The data file has to have a symbol column for this to work correctly.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT combined.tsv
# PARAMETER produce.identifiers [no, yes] DEFAULT no (Should approximate identifiers for the gene be returned)


# Combine probes or probeset to genes
# 18.12.2008 JTT
#
# modified 3.3.2010, MG, to include symbol and description columns in the output
# if present in the input data table

# Loads the file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Check whether there are symbol and description columns
symbols <- dat[,grep("symbol", names(dat))]
descriptions <- dat[,grep("description", names(dat))]

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
	
	# Fetch back the symbol and description column, if present
	# in the original data
	if (length(symbols)>0) {
		symbols2 <- as.character(dat$symbol[match (rownames(dat6), rownames(dat))])
		descriptions2 <- as.character(dat$description[match (rownames(dat6), rownames(dat))])
		annotations <- cbind(symbols2, descriptions2)
		dat6 <- data.frame(annotations, dat6)
		names(dat6) [1:2] <- c("symbol", "description")
	}
	
	write.table(data.frame(dat6), file="combined.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	
} else {
	
	# Putting a data file together
	colnames(dat3)<-colnames(dat2)
	# colnames(dat4)<-paste("sd.", colnames(dat2), sep="")
	symbol<-m$Group.1
	# dat5<-data.frame(symbol, dat3, dat4) 
	dat5<-data.frame(dat3) 
	rownames(dat5)<-symbol
	
	# Fetch back the description column, if present
	# in the original data
	if (length(descriptions)>0) {
		descriptions2 <- as.character(dat$description[match (rownames(dat5), dat$symbol)])
		dat5 <- data.frame(descriptions2, dat5)
		names(dat5) [1] <- "description"
	}
	
	# Writing data to disk
	write.table(data.frame(dat5), file="combined.tsv", col.names=T, quote=F, sep="\t", row.names=T)
}
