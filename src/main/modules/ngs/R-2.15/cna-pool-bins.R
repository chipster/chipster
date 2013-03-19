# TOOL cna-pool-bins.R: "Pool binned read counts" (Pools binned read counts together from two or more files.)
# INPUT binned-hits-{...}.tsv: "Individual files with binned hits" TYPE GENERIC
# OUTPUT binned-hits-pooled.tsv: "Pooled read counts"

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-06-07

filenames <- list.files(pattern='^binned-hits-[0-9]*\\.tsv$')

dat <- read.table(filenames[1], header=TRUE, sep='\t', row.names=1, as.is=TRUE)

if (length(filenames) > 1) {
  for (i in 2:length(filenames)) {
    x <- read.table(filenames[i], header=TRUE, sep='\t', row.names=1, as.is=TRUE)
    if (nrow(x) != nrow(dat))
      stop("CHIPSTER-NOTE: All input files need to be binned using the same bin size.")
    dat$count <- dat$count + x$count
    dat$corrected <- dat$corrected + x$corrected
  }
}

options(scipen=10)
write.table(dat, 'binned-hits-pooled.tsv', quote=FALSE, sep='\t', na='')

# EOF
