# TOOL cna-define-experiment.R: "Define CNA-seq experiment" (This tool creates a phenodata file containing descriptive information about samples and experiment setup.)
# INPUT binned-hits-{...}.tsv: "Individual files with binned hits" TYPE GENERIC
# OUTPUT cna-data-table.tsv: "Data table with log-transformed read counts"
# OUTPUT META phenodata.tsv: "Experiment description file"
# PARAMETER counts: "Counts" TYPE [count: "Original raw counts", corrected: "GC corrected counts"] DEFAULT corrected (Whether to use original raw read counts, or GC corrected ones.)
# PARAMETER normalization: "Normalization" TYPE [median: median] DEFAULT median (Normalization method.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-12-22

source(file.path(chipster.common.path, 'MPScall', 'CGHcallPlus-R-2.12.R'))

if (counts == 'count') {
  colClasses1 <- c('character', 'character', 'integer', 'integer', 'numeric', 'NULL')
  colClasses2 <- c('character', 'NULL', 'NULL', 'NULL', 'numeric', 'NULL')
} else {
  colClasses1 <- c('character', 'character', 'integer', 'integer', 'NULL', 'numeric')
  colClasses2 <- c('character', 'NULL', 'NULL', 'NULL', 'NULL', 'numeric')
}

filenames <- list.files(pattern='^binned-hits-[0-9]*\\.tsv$')

dat <- read.table(filenames[1], header=TRUE, sep='\t', row.names=1, colClasses=colClasses1)
dat[,4] <- dat[,4] - min(dat[,4], na.rm=TRUE) + 1 # to prevent negative values
dat[,4] <- log2(dat[,4])
dat <- data.frame(bin=1:nrow(dat), dat)

if (length(filenames) > 1) {
  for (i in 2:length(filenames)) {
    x <- read.table(filenames[i], header=TRUE, sep='\t', row.names=1, colClasses=colClasses2)
    if (nrow(x) != nrow(dat))
      stop("CHIPSTER-NOTE: All input files need to be binned using the same bin size.")
    x[,1] <- x[,1] - min(x[,1], na.rm=TRUE) + 1 # to prevent negative values
    dat <- cbind(dat, log2(x[,1]))
  }
}

dat$chromosome[dat$chromosome=='X'] <- '23'
dat$chromosome[dat$chromosome=='Y'] <- '24'
dat$chromosome[dat$chromosome=='MT'] <- '25'
dat$chromosome <- as.integer(dat$chromosome)

cgh.raw <- make_cghRaw(dat)
cgh.pre <- preprocess(cgh.raw, nchrom=max(chromosomes(cgh.raw)))
cgh.nor <- normalize(cgh.pre, method=normalization)

dat2 <- data.frame(cgh.nor@featureData@data, round(copynumber(cgh.nor), digits=2))
colnames(dat2) <- c('chromosome', 'start', 'end', sprintf('chip.microarray%.3i', 1:length(filenames)))

dat2$chromosome <- as.character(dat2$chromosome)
dat2$chromosome[dat2$chromosome=='23'] <- 'X'
dat2$chromosome[dat2$chromosome=='24'] <- 'Y'
dat2$chromosome[dat2$chromosome=='25'] <- 'MT'

# generate phenodata

phenodata <- data.frame(sample=filenames, chiptype='not applicable', experiment='cna_seq', group='', stringsAsFactors=FALSE) # description?

# write outputs

options(scipen=10)
write.table(dat2, 'cna-data-table.tsv', quote=FALSE, sep='\t', na='')

write.table(phenodata, file='phenodata.tsv', quote=FALSE, sep='\t', row.names=FALSE)

# EOF
