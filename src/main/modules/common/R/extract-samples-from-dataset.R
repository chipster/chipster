# TOOL extract-samples-from-dataset.R: "Extract samples from dataset" (Makes a subset of a dataset. The samples to be extracted have to be coded with 1, and the samples to be deleted with 0 in the same phenodata column. If there are missing values in the specified phenodata column, the samples that do have a value are extracted.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT extract.tsv: extract.tsv 
# OUTPUT META phenodata-extract.tsv: phenodata-extract.tsv 
# PARAMETER column.extract: "Which column contains the extraction info"  TYPE METACOLUMN_SEL DEFAULT group (Phenodata column containing the samples to be extracted)

# JTT, 19.10.2007
# MG 20.4.2010, ability to include annotation info in the output dataset  
# IS 28.7.2010, modified to first filter out samples with missing values. And in the case of no missing values, follow the earlier behaviour of coding samples to be extracted with 1 and samples to be removed with 0. 
# MG 17.6.2011, ability to handle 2-color array data when extracting single sample. 
# IS 18.3.2013, modified to be able to handle any matrices (chip, flag, ...) and annotation columns present in the data.
# MK 30.4.2013, added the ability to deal with NGS count tables
# IS 6.1.2014, making less sensitive about phenodata columnn names

# Loads the data file
file <- c('normalized.tsv')
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, as.is=TRUE, check.names=FALSE)

# Loads phenodata
phenodata <- read.table('phenodata.tsv', header=TRUE, sep='\t', as.is=TRUE, check.names=FALSE)

# Extract the data from the phenodata column
extract <- as.vector(phenodata[,column.extract])
extract[extract == ''] <- NA

# If there are samples with missing values, extract the ones that do have values.
if (sum(is.na(extract)) > 0) {
  extract <- ifelse(is.na(extract), 0, 1)
}

extract <- as.integer(extract)

# Sanity checks
if(length(unique(extract))>2) {
	stop("CHIPSTER-NOTE: You have specified more than two groups! You need to define exactly two groups.")
}
if(max(extract)>1) {
	stop("CHIPSTER-NOTE: The groups should be defined with values of 0 and 1! You have numbers larger than 1 in the definitions.")
}

# identify different matrices (chip, flag, ...) present in the data
x <- colnames(dat)
suffix <- sub('^chip\\.', '', x[grep('^chip\\.', x)[1]]) #pick up the name of the first experiment
matrices <- sub(suffix, '', x[grep(paste("\\.", suffix, sep=""), x)]) #remove the suffix and find all columns having the name
annotations <- 1:length(x)
for (m in matrices) {
	annotations <- setdiff(annotations, grep(m, x))
}


if(length(annotations) > 0) {
	dat2 <- dat[,annotations]
} else {
	dat2 <- data.frame(row.names=rownames(dat)) 
}
for (m in matrices) {
	dat2 <- cbind(dat2, dat[,grep(m, x), drop=FALSE][,extract==1, drop=FALSE])
}


phenodata2<-phenodata[which(extract==1),]

# update aberration frequencies if present
if ('loss.freq' %in% x) {
  calls <- dat2[,grep('^flag\\.', colnames(dat2)), drop=FALSE]
  dat2$loss.freq <- round(rowMeans(calls < 0), digits=3)
  dat2$gain.freq <- round(rowMeans(calls == 1), digits=3)
  if (2 %in% calls) {
    dat2$amp.freq <- round(rowMeans(calls == 2), digits=3)
  } else {
    dat2$amp.freq <- NULL
  }
}

# Writing the data to disk
options(scipen=10)
write.table(dat2, file="extract.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(phenodata2, file="phenodata-extract.tsv", sep="\t", row.names=F, col.names=T, quote=F, na='')


# EOF
