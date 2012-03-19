# TOOL cna-define-experiment.R: "Define CNA-seq experiment" (This tool creates a phenodata file containing descriptive information about samples and experiment setup.)
# INPUT binned-hits-{...}.tsv: "Individual files with binned hits" TYPE GENERIC
# OUTPUT cna-data-table.tsv: "Data table with log-transformed read counts"
# OUTPUT META phenodata.tsv: "Experiment description file"
# PARAMETER counts: "Counts" TYPE [count: "Original raw counts", corrected: "GC corrected counts"] DEFAULT count (Whether to use original raw read counts, or GC corrected ones.)
# PARAMETER log2: "Log2 transform counts" TYPE [no: no, yes: yes] DEFAULT no (Whether the counts should be log2 transformed.)
# PARAMETER normalization: "Normalization" TYPE [none: none, median: median] DEFAULT none (Normalization method.)
# PARAMETER min.mappability: "Mimimum mappability" TYPE DECIMAL FROM 0 TO 100 DEFAULT 0 (The bins with lower mappability will be removed.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-03-01

source(file.path(chipster.tools.path, 'MPScall', 'CGHcallPlus-R-2.12.R'))

colClasses1 <- c('character', 'character', 'integer', 'integer', 'numeric', 'numeric')
colClasses2 <- c('character', 'NULL', 'NULL', 'NULL', 'numeric', 'numeric')

filenames <- list.files(pattern='^binned-hits-[0-9]*\\.tsv$')

dat <- read.table(filenames[1], header=TRUE, sep='\t', row.names=1, colClasses=colClasses1)
total.reads <- sum(dat[,4], na.rm=TRUE)
if (count == 'corrected') {
  dat[,4] <- NULL
} else {
  dat[,5] <- NULL
}
if (log2 == 'yes')
  dat[,4] <- log2(dat[,4] - min(dat[,4], na.rm=TRUE) + 1) # to prevent negative values

if (length(filenames) > 1) {
  for (i in 2:length(filenames)) {
    x <- read.table(filenames[i], header=TRUE, sep='\t', row.names=1, colClasses=colClasses2)
    if (nrow(x) != nrow(dat))
      stop("CHIPSTER-NOTE: All input files need to be binned using the same bin size.")
    total.reads <- c(total.reads, sum(dat[,1], na.rm=TRUE))
    if (count == 'corrected') {
      dat[,1] <- NULL
    } else {
      dat[,2] <- NULL
    }
    if (log2 == 'yes')
      x[,1] <- log2(x[,1] - min(x[,1], na.rm=TRUE) + 1) # to prevent negative values
    dat <- cbind(dat, x[,1])
  }
}

identifiers <- gsub('tsv$', '', filenames)
identifiers <- gsub('\\.', '', identifiers)
identifiers <- gsub('-', '', identifiers)
colnames(dat)[-(1:3)] <- paste('chip.', identifiers, sep='')

if (min.mappability > 0) {
  mappability <- read.table(file.path(chipster.tools.path, 'MPScall', genome.build, paste('mappability.', bin.size, 'kbp.txt.gz', sep='')), header=TRUE, sep='\t', as.is=TRUE, colClasses=c('character', 'integer', 'integer', 'numeric'))
  dat <- dat[mappability$mappability >= min.mappability,]
}

if (normalization != 'none') {
  dat$chromosome[dat$chromosome=='X'] <- '23'
  dat$chromosome[dat$chromosome=='Y'] <- '24'
  dat$chromosome[dat$chromosome=='MT'] <- '25'
  dat$chromosome <- as.integer(dat$chromosome)

  cgh.raw <- make_cghRaw(dat)
  cgh.pre <- preprocess(cgh.raw, nchrom=max(chromosomes(cgh.raw)))
  cgh.nor <- normalize(cgh.pre, method=normalization)

  dat2 <- data.frame(cgh.nor@featureData@data, round(copynumber(cgh.nor), digits=2))
  colnames(dat2) <- c('chromosome', 'start', 'end', paste('chip.', identifiers, sep=''))

  dat2$chromosome <- as.character(dat2$chromosome)
  dat2$chromosome[dat2$chromosome=='23'] <- 'X'
  dat2$chromosome[dat2$chromosome=='24'] <- 'Y'
  dat2$chromosome[dat2$chromosome=='25'] <- 'MT'
  dat <- dat2
}

# generate phenodata

phenodata <- data.frame(sample=filenames, chiptype='not applicable', experiment='cna_seq', reads=total.reads, group='', stringsAsFactors=FALSE)

# write outputs

options(scipen=10)
write.table(dat, 'cna-data-table.tsv', quote=FALSE, sep='\t', na='')

write.table(phenodata, file='phenodata.tsv', quote=FALSE, sep='\t', row.names=FALSE)

# EOF
