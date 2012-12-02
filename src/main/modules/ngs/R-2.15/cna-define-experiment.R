# TOOL cna-define-experiment.R: "Define CNA-seq experiment" (This tool creates a phenodata file containing descriptive information about samples and experiment setup.)
# INPUT binned-hits-{...}.tsv: "Individual files with binned hits" TYPE GENERIC
# OUTPUT cna-data-table.tsv: "Data table with log-transformed read counts"
# OUTPUT META phenodata.tsv: "Experiment description file"
# PARAMETER counts: "Read counts" TYPE [count: "Original raw counts", corrected: "GC corrected counts"] DEFAULT count (Whether to use original raw read counts, or GC corrected ones.)
# PARAMETER log2.transformation: "log2 transform counts" TYPE [no: no, yes: yes] DEFAULT no (Whether the counts should be log2 transformed.)
# PARAMETER min.mappability: "mimimum mappability" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0 (The bins with lower mappability will be removed. Values are between 0 and 1.)
# PARAMETER normalization: "normalization" TYPE [none: none, median: median] DEFAULT none (Normalization method.)
# PARAMETER organism: "organism" TYPE [human: human] DEFAULT human (Organism.)
# PARAMETER genome.build: "human genome build" TYPE [GRCh37: GRCh37] DEFAULT GRCh37 (Genome build.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-03-20

source(file.path(chipster.common.path, 'CGHcallPlus.R'))

filenames <- list.files(pattern='^binned-hits-[0-9]*\\.tsv$')
dat <- read.table(filenames[1], header=TRUE, sep='\t', row.names=1, as.is=TRUE)
total.reads <- sum(dat$count, na.rm=TRUE)
if (counts == 'corrected') {
  dat$count <- NULL
} else {
  dat$corrected <- NULL
}

if (length(filenames) > 1) {
  for (i in 2:length(filenames)) {
    x <- read.table(filenames[i], header=TRUE, sep='\t', row.names=1, as.is=TRUE)
    if (nrow(x) != nrow(dat))
      stop("CHIPSTER-NOTE: All input files need to be binned using the same bin size.")
    total.reads <- c(total.reads, sum(x$count, na.rm=TRUE))
    dat <- cbind(dat, x[,counts])
  }
}

identifiers <- gsub('tsv$', '', filenames)
identifiers <- gsub('\\.', '', identifiers)
identifiers <- gsub('-', '', identifiers)
colnames(dat)[-(1:3)] <- paste('chip.', identifiers, sep='')

if (log2.transformation == 'yes') {
  for (i in colnames(dat)[-(1:3)])
    dat[,i] <- log2(dat[,i] - min(dat[,i], na.rm=TRUE) + 1)
}

if (min.mappability > 0) {
  bin.size <- (dat$end[1] - dat$start[1] + 1) / 1000
  min.mappability <- 100 * min.mappability
  mappability <- read.table(file.path(chipster.tools.path, 'MPScall', genome.build, paste('mappability.', bin.size, 'kbp.txt.gz', sep='')), header=TRUE, sep='\t', as.is=TRUE, colClasses=c('character', 'integer', 'integer', 'numeric'))
  dat <- dat[mappability$mappability >= min.mappability,]
}

if (normalization != 'none') {
  dat2 <- data.frame(bin=rownames(dat), dat, stringsAsFactors=FALSE)
  colnames(dat2) <- c('bin', colnames(dat))
  dat2$chromosome[dat2$chromosome=='X'] <- '23'
  dat2$chromosome[dat2$chromosome=='Y'] <- '24'
  dat2$chromosome[dat2$chromosome=='MT'] <- '25'
  dat2$chromosome <- as.integer(dat2$chromosome)

  cgh.raw <- make_cghRaw(dat2)
  cgh.pre <- preprocess(cgh.raw, nchrom=24)
  cgh.nor <- normalize(cgh.pre, method=normalization)

  dat <- data.frame(cgh.nor@featureData@data, round(copynumber(cgh.nor), digits=2))
  colnames(dat) <- c('chromosome', 'start', 'end', paste('chip.', identifiers, sep=''))

  dat$chromosome <- as.character(dat$chromosome)
  dat$chromosome[dat$chromosome=='23'] <- 'X'
  dat$chromosome[dat$chromosome=='24'] <- 'Y'
  dat$chromosome[dat$chromosome=='25'] <- 'MT'
}

# generate phenodata

phenodata <- data.frame(sample=filenames, chiptype='not applicable', experiment='cna_seq', reads=total.reads, group='', stringsAsFactors=FALSE)

# write outputs

options(scipen=10)
write.table(dat, 'cna-data-table.tsv', quote=FALSE, sep='\t', na='')

write.table(phenodata, file='phenodata.tsv', quote=FALSE, sep='\t', row.names=FALSE)

# EOF
