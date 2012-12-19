# TOOL cna-define-experiment.R: "Define CNA-seq experiment" (This tool creates a phenodata file containing descriptive information about samples and experiment setup.)
# INPUT binned-hits-{...}.tsv: "Individual files with binned hits" TYPE GENERIC
# OUTPUT cna-data-table.tsv: "Data table with log-transformed read counts"
# OUTPUT META phenodata.tsv: "Experiment description file"

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-12-18

# look for sample names? !!!

filenames <- list.files(pattern='^binned-hits-[0-9]*\\.tsv$')
dat <- read.table(filenames[1], header=TRUE, sep='\t', row.names=1, as.is=TRUE)
total.reads <- sum(dat$count, na.rm=TRUE)
dat$corrected <- NULL
dat$count <- log2(dat$count+1)

if (length(filenames) > 1) {
  for (i in 2:length(filenames)) {
    x <- read.table(filenames[i], header=TRUE, sep='\t', row.names=1, as.is=TRUE)
    if (nrow(x) != nrow(dat))
      stop("CHIPSTER-NOTE: All input files need to be binned using the same bin size.")
    total.reads <- c(total.reads, sum(x$count, na.rm=TRUE))
    x$corrected <- NULL
    dat <- cbind(dat, log2(x$count+1))
  }
}

identifiers <- gsub('tsv$', '', filenames)
identifiers <- gsub('\\.', '', identifiers)
identifiers <- gsub('-', '', identifiers)
colnames(dat)[-(1:3)] <- paste('chip.', identifiers, sep='')

# generate phenodata
phenodata <- data.frame(sample=filenames, chiptype='not applicable', experiment='cna_seq', reads=total.reads, group='', stringsAsFactors=FALSE)

# write outputs
options(scipen=10)
write.table(dat, 'cna-data-table.tsv', quote=FALSE, sep='\t', na='')
write.table(phenodata, file='phenodata.tsv', quote=FALSE, sep='\t', row.names=FALSE)

# EOF
