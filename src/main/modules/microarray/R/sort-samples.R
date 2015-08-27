# TOOL sort-samples.R: "Sort samples" (Sorts samples according to a phenodata column. The column should contain numerical values, since the samples are sorted in ascending order according to the values.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT sort-samples.tsv: sort-samples.tsv 
# OUTPUT META phenodata-sorted.tsv: phenodata-sorted.tsv 
# PARAMETER column: Column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column specifying how to sort)

# JTT, 06.02.2008, Sort samples
# MG,  16.11.2010, modified to also generate a re-ordered phenodata file to reflect the re-ordered data

# Default parameters
#column<-"group"
 
# Loads libraries
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Loads data (which file to search)
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Extracts the column from phenodata
groups<-phenodata[,pmatch(column,colnames(phenodata))]

# identify matrices (chip, flag, segmented, ...) present in the data
datnames <- colnames(dat)
suffix <- sub('^chip\\.', '', datnames[grep('^chip\\.', datnames)[1]])	#name of the first array with the chip-prefix
matrices <- sub(suffix, '', datnames[grep(suffix, datnames)])			#returns chip.,  flags., etc

#create data matrices
dat3 <- dat
for (m in matrices) {
	dat_temp <- dat[, grep(m, datnames)]
	dat_temp <- dat_temp[,order(groups)]	
	dat3[, grep(m, datnames)] <- dat_temp
	names(dat3)[grep(m, datnames)] <- names(dat_temp)
}

phenodata2<-phenodata[order(groups),]

write.table(dat3, file="sort-samples.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(phenodata2, file="phenodata-sorted.tsv", sep="\t", row.names=F, col.names=T, quote=F)
