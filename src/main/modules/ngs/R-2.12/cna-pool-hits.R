# TOOL cna-compare-to-reference.R: "Compare to reference sample" (This tool compares binned read counts of a test sample to binned read counts of a reference sample by taking a ratio between the two.)
# INPUT binned-hits-{...}.tsv: "Individual files with binned hits" TYPE GENERIC
# OUTPUT binned-hits-pooled.tsv: "Pooled read counts"

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-02-15

filenames <- list.files(pattern='^binned-hits-[0-9]*\\.tsv$')

dat <- read.table(filenames[1], header=TRUE, sep='\t', row.names=1, colClasses=c('character', 'character', 'integer', 'integer', 'numeric', 'numeric'))
dat[,4] <- dat[,4] - min(dat[,4], na.rm=TRUE) + 1 # to prevent negative values
dat[,4] <- log2(dat[,4])
dat <- data.frame(bin=1:nrow(dat), dat)

if (length(filenames) > 1) {
  for (i in 2:length(filenames)) {
    x <- read.table(filenames[i], header=TRUE, sep='\t', row.names=1, colClasses=c('character', 'NULL', 'NULL', 'NULL', 'numeric', 'numeric'))
    if (nrow(x) != nrow(dat))
      stop("CHIPSTER-NOTE: All input files need to be binned using the same bin size.")
    dat[,4] <- dat[,4] + x[,1]
    dat[,5] <- dat[,5] + x[,2]
  }
}

options(scipen=10)
write.table(dat, 'binned-hits-pooled.tsv', quote=FALSE, sep='\t', na='')

# EOF