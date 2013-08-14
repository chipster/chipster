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

# Identify replicates
del.cols   <- phenodata[, column1]
del.groups <- unique(phenodata[which(del.cols == 1), column2])
sample.groups <- phenodata[, column2]

dat.sub <- matrix(0, ncol=ncol(dat), nrow(dat));
for (i in 1:length(del.groups)) {
	if(length(which(del.cols %in% 1 & sample.groups %in% del.groups[i] == TRUE)) == 1) {
		dat.sub[, which(sample.groups == del.groups[i])] <- dat[, (del.cols %in% 1 & sample.groups %in% del.groups[i])]
	} else {
		dat.sub[, which(sample.groups == del.groups[i])] <- apply(dat[, (del.cols %in% 1 & sample.groups %in% del.groups[i])], 1, mean)
	}

}

# 
dat2 <- dat - dat.sub
dat2 <- dat2[, -(which(del.cols==1))]

phenodata2 <- phenodata[-(which(del.cols==1)),]
phenodata2$sample <- sprintf('microarray%.3i', 1:nrow(phenodata2))
colnames(dat2) <- phenodata2$sample
#rownames(phenodata2) <- sprintf('microarray%.3i', 1:nrow(phenodata2))

# Writes out the combined table
write.table(dat2, "deleted.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(phenodata2, file='phenodata.tsv', quote=FALSE, sep='\t', row.names=FALSE)

c