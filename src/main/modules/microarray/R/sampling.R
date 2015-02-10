# TOOL sampling.R: "Random sampling" (Generates a random sample of the data.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT sample.tsv: sample.tsv 
# PARAMETER what.to.sample: "What to sample" TYPE [genes: genes, chips: chips] DEFAULT genes (Sample genes or chips)
# PARAMETER sample.size: "Sample size" TYPE DECIMAL FROM 0 TO 100 DEFAULT 10 (Percentage of the data to sample)
# PARAMETER replacement: Replacement TYPE [yes: yes, no: no] DEFAULT no (Sampling with or without replacement)

# JTT: 20.10.2008, Random sampling
# MK: 31.10.2013, sampling with replacement option added

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Sanity checks
if(ncol(dat)<2 | nrow(dat)<2) {
   stop("CHIPSTER-NOTE: You need to have at least 2 data points from which to sample!")
}

# identify matrices (chip, flag, segmented, ...) present in the data
datnames <- colnames(dat)
suffix <- sub('^chip\\.', '', datnames[grep('^chip\\.', datnames)[1]])	#name of the first array with the chip-prefix
matrices <- sub(suffix, '', datnames[grep(suffix, datnames)])			#returns chip.,  flags., etc

if(replacement == "yes") {
	replacement <- TRUE
} else {
	replacement <- FALSE
}

# Sampling
if(what.to.sample=="genes") {
   v <- sample(1:nrow(dat), floor(nrow(dat)*sample.size/100), replace=replacement)
   dat2<-dat[v,]
}

if(what.to.sample=="chips") {
    v <- sample(1:length(grep(matrices[1], datnames)) , floor(length(grep(matrices[1], datnames))*sample.size/100), replace=replacement)

	annotations <- 1:ncol(dat)
	for (m in matrices) {
		annotations <- setdiff(annotations, grep(m, datnames))
	}
	dat2 <- dat[,annotations]

	for (m in matrices) {
		dat_temp <- dat[, grep(m, datnames)]
		dat_temp <- dat_temp[, v]
		dat2 <- cbind(dat2, dat_temp)
	}
}

# Writing data to disk
write.table(data.frame(dat2), file="sample.tsv", col.names=T, quote=F, sep="\t", row.names=T)
