# ANALYSIS "aCGH tools (beta testing)"/"Count overlapping CNVs" (Counts overlapping CNVs from the database of genomic variants.)
# INPUT GENERIC normalized.tsv
# OUTPUT normalized.tsv
# PARAMETER genome.build [GRCh37, NCBI36, NCBI35] DEFAULT GRCh37 (The genome build to use for fetching the CNV data.)

# count-overlapping-cnvs.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-08-13

dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

first.data.col <- min(grep('chip', names(dat)), grep('flag', names(dat)))

dat2 <- dat[1:first.data.col-1]

# load cnvs
if (genome.build=='NCBI35') {
  cnv <- read.table('http://projects.tcag.ca/variation/downloads/variation.hg17.v9.mar.2010.txt', header=TRUE, sep='\t', as.is=TRUE)
} else if (genome.build=='NCBI36') {
  cnv <- read.table('http://projects.tcag.ca/variation/downloads/variation.hg18.v9.mar.2010.txt', header=TRUE, sep='\t', as.is=TRUE)
} else {
  cnv <- read.table('http://projects.tcag.ca/variation/downloads/variation.hg19.v9.mar.2010.txt', header=TRUE, sep='\t', as.is=TRUE)
}
cnv <- cnv[cnv$VariationType == 'CopyNumber',]
cnv$Chr <- substr(cnv$Chr, 4, 5)

dat2$cnvs <- 0
for (i in rownames(dat)) {
  dat2[i, 'cnvs'] <- nrow(subset(cnv, Chr == dat[i, 'chromosome'] & 
                                    Start <= dat[i, 'end'] & 
                                      End >= dat[i, 'start']))
}

# calculate density of CNVs per Mb
dat2$cnv.density <- dat2$cnvs / (dat2$end - dat2$start + 1) * 1000000

dat2 <- cbind(dat2, dat[,first.data.col:ncol(dat)])

write.table(dat2, file='normalized.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF