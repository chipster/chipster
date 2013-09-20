# TOOL merge-pheno-and-data.R: "Merge pheno- and expression matrices" (Merge phenodata and expression matrices into a single spreadsheet for easy data import.)
# INPUT normalized.tsv TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT merged.tsv
# PARAMETER excel.file TYPE [yes, no] DEFAULT no (Fix column headers so that the file opens correctly in Excel)

# MK 20.09.2013

# Sanity checks
file<-c("phenodata.tsv")
dat1<-read.table(file, sep="\t", header=T, row.names=1,  quote = "")

file<-c("normalized.tsv")
dat2<-read.table(file, sep="\t", header=T, row.names=1,  quote = "")

#create an empty matrix
temp <- data.frame(matrix(NA, ncol=ncol(dat2), nrow=ncol(dat1)+1))
colnames(temp) <- colnames(dat2)
dat3 <- rbind(temp, dat2)

#find expression columns associated with the first sample and check what it their prefix
firstsample <- colnames(dat2)[grep(rownames(dat1)[1], colnames(dat2))]
prefixes 	<- gsub(rownames(dat1)[1], "", firstsample)

for(i in 1:length(prefixes)) {
	temp <- dat1[!is.na(match(paste(prefixes[i],rownames(dat1),sep=""), colnames(dat2))), ]
	cols <- match(paste(prefixes[i],rownames(temp),sep=""), colnames(dat2))
	dat3[1, cols] <- rownames(temp)
	dat3[2:(ncol(dat1)+1), cols] <- t(temp)
}

rownames(dat3) <- make.names(c("sample", colnames(dat1), rownames(dat2)), unique=TRUE)

# Writes out the combined table
if(excel.file == "no") {
	write.table(dat3, "merged.tsv", sep="\t", row.names=T, col.names=T, quote=F)
} else {
	write.table(dat3, "merged.tsv", sep="\t", row.names=T, col.names=NA, quote=F)
}