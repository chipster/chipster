# ANALYSIS "aCGH tools (beta testing)"/"Merge data sets" (Merges two data sets. Only rows common to both are kept. All annotation data is taken from the first data set.)
# INPUT GENE_EXPRS normalized_1.tsv, GENE_EXPRS normalized_2.tsv, GENERIC phenodata_1.tsv, GENERIC phenodata_2.tsv
# OUTPUT merged.tsv, phenodata-merged.tsv

# merge-datasets.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-08-19

# read data set 1
dat <- read.table('normalized_1.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

# read data set 2
dat2 <- read.table('normalized_2.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

# take common rows
common.rows <- intersect(rownames(dat), rownames(dat2))
dat <- dat[common.rows,]
dat2 <- dat2[common.rows,]

# extract data
ratios <- dat[,grep('chip', colnames(dat))]
calls <- dat[,grep('flag', colnames(dat))]
ratios <- cbind(ratios, dat2[,grep('chip', colnames(dat2))])
calls <- cbind(calls, dat2[,grep('flag', colnames(dat2))])

# generate new identifiers
microarrays <- sprintf('microarray%.3i', 1:ncol(ratios))
colnames(ratios) <- paste('chip.', microarrays, sep='')
if (ncol(ratios) == ncol(calls))
  colnames(calls) <- paste('flag.', microarrays, sep='')

# remove old matrices from data set 1
dat <- dat[,-grep('chip', colnames(dat))]
if (length(grep('flag', colnames(dat))) > 0)
  dat <- dat[,-grep('flag', colnames(dat))]

# calculate new frequencies
if ('loss.freq' %in% colnames(dat) && ncol(ratios) == ncol(calls)) {
  dat$loss.freq <- round(mean(as.data.frame(t(calls==-1))), digits=3)
  dat$gain.freq <- round(mean(as.data.frame(t(calls==1))), digits=3)
  if (2 %in% calls) {
    dat$amp.freq <- round(mean(as.data.frame(t(calls==2))), digits=3)
  } else {
    dat$amp.freq <- NULL
  }
} else {
  dat$loss.freq <- NULL
  dat$gain.freq <- NULL
  dat$amp.freq <- NULL
}

# generate new table
dat <- cbind(dat, ratios)
if (ncol(ratios) == ncol(calls))
  dat <- cbind(dat, calls)

# process phenodata
phenodata1 <- read.table('phenodata_1.tsv', header=TRUE, sep='\t', as.is=TRUE)
phenodata2 <- read.table('phenodata_2.tsv', header=TRUE, sep='\t', as.is=TRUE)

# fill in columns present only in one phenodata table
for (col in setdiff(colnames(phenodata1), colnames(phenodata2)))
  phenodata2[,col] <- NA
for (col in setdiff(colnames(phenodata2), colnames(phenodata1)))
  phenodata1[,col] <- NA

# combine phenodata tables and update sample identifiers
phenodata <- rbind(phenodata1, phenodata2)
phenodata$sample <- microarrays

# write files
write.table(dat, file='merged.tsv', quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
write.table(phenodata, file='phenodata-merged.tsv', quote=FALSE, sep='\t', na='', row.names=FALSE, col.names=TRUE)

# EOF