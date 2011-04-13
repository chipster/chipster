# TOOL search-coordinates.R: "Search by genomic coordinates" (Search by chromosome name and starting and ending base pair positions. The data must contain the corresponding columns.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT search-coordinates.tsv: search-coordinates.tsv 
# PARAMETER position: position STRING (Position to search for. Must contain three values that are separated by tabs, hyphens, colons or two dots (e.g. X:100-200 or 7:600..700\).)
# PARAMETER include.partial.overlaps: include.partial.overlaps [yes: yes, no: no] DEFAULT yes (Whether to include only features that are completely contained within the search window, or also partial overlaps.)

# search-coordinates.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-04-13

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
write.table(format(dat, scientific=FALSE), file='search-coordinates.tsv', quote=FALSE, sep='\t')

# EOF
