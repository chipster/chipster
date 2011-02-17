# ANALYSIS Utilities/"Search by genomic coordinates" ()
# INPUT GENE_EXPRS normalized.tsv
# OUTPUT search-coordinates.tsv
# PARAMETER chromosome STRING DEFAULT 1 ()
# PARAMETER start INTEGER DEFAULT 100000000 ()
# PARAMETER end INTEGER DEFAULT 0 ()
# PARAMETER number.of.closest.results.to.return INTEGER DEFAULT 5 ()

# search-coordinates.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-11-17

# load inputs
dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

# discard data from other chromosomes
dat <- dat[dat$chromosome == chromosome,]

# add distance column
dat2 <- dat[,pos]
dat2$distance <- abs(dat$start + (dat$end - dat$start) / 2 - (start + (end - start) / 2))
dat2 <- cbind(dat2, dat[,setdiff(colnames(dat), pos)])
dat <- dat2

# order according to distance
dat <- dat[order(dat$distance),]

# check that we are not trying to return more results than we have
number.of.closest.results.to.return <- min(number.of.closest.results.to.return, nrow(dat))

if (number.of.closest.results.to.return == 0) {
  dat <- dat[0,]
} else {
  dat <- dat[1:number.of.closest.results.to.return,]
}

# write output
write.table(dat, file='search-coordinates.tsv', quote=FALSE, sep='\t')

# EOF
