# TOOL delete-and-subtract-columns.R: "Delete and subtract columns" (Delete the specified column or columns from the data and subtract their values from associated samples.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT deleted.tsv: deleted.tsv
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER column1: Delete TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing the samples to be deleted. Chips to be deleted should be coded with 1 in the column.)
# PARAMETER column2: "Sample Groups" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing sample groups. Effects of deleted chips of one group will be removed from the effects of the other chips of that group. If more than one chip are to be removed from a group, removed chips will be averaged.)

# Deletes columns from data and subtacts their (average values) from other columns which they are associated with
# MK 13.08.2013

# Sanity checks
if(column1=="EMPTY" || column2=="EMPTY") {
	stop("You need to select columns describing the samples to delete and sample groups!")
}

# Loads the data
file<-c("normalized.tsv")
dat<-read.table(file, sep="\t", header=T, row.names=1)
phenodata <- read.table('phenodata.tsv', header=TRUE, sep='\t', as.is=TRUE)

# identify matrices (chip, flag, segmented, ...) present in the data
datnames <- colnames(dat)
suffix <- sub('^chip\\.', '', datnames[grep('^chip\\.', datnames)[1]])	#name of the first array with the chip-prefix
matrices <- sub(suffix, '', datnames[grep(suffix, datnames)])			#returns chip.,  flags., etc

# identify annotation columns (that are not part of any of the matrices)
annotations <- 1:ncol(dat)
for (m in matrices) {
	annotations <- setdiff(annotations, grep(m, datnames))
}
dat2 <- dat[,annotations]

#identify samples to be deleted
del.cols   <- phenodata[, column1]
del.groups <- unique(phenodata[which(del.cols == 1), column2])
sample.groups <- phenodata[, column2]

#create new phenodata
phenodata2 <- phenodata[-(which(del.cols==1)),]
phenodata2$sample <- sprintf('microarray%.3i', 1:nrow(phenodata2))
rownames(phenodata2) <- sprintf('microarray%.3i', 1:nrow(phenodata2))

#create data matrices
for (m in matrices) {
	dat3 <- dat[, grep(m, datnames)]
	if(m == "chip.") {
		dat.sub <- matrix(0, ncol=ncol(dat3), nrow(dat3));
		for (i in 1:length(del.groups)) {
			if(length(which(del.cols %in% 1 & sample.groups %in% del.groups[i] == TRUE)) == 1) {
				dat.sub[, which(sample.groups == del.groups[i])] <- dat3[, (del.cols %in% 1 & sample.groups %in% del.groups[i])]
			} else {
				dat.sub[, which(sample.groups == del.groups[i])] <- apply(da3t[, (del.cols %in% 1 & sample.groups %in% del.groups[i])], 1, mean)
			}
		}
		dat4 <- dat3 - dat.sub
		dat4 <- dat4[, -(which(del.cols==1))]
	} else {		
		#handle flag information
		dat4 <- dat3[, -(which(del.cols==1))]
	}
	colnames(dat4) <- paste(m, phenodata2$sample, sep="")
	dat2 <- cbind(dat2, dat4)
}		

# Writes out the combined table
write.table(dat2, file="deleted.tsv", sep="\t", quote=F)
write.table(phenodata2, file='phenodata.tsv', quote=FALSE, sep='\t', row.names=FALSE)
