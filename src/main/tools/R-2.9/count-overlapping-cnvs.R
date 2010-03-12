# ANALYSIS "aCGH tools (beta testing)"/"Count overlapping CNVs" (Counts overlapping CNVs from the database of genomic variants.)
# INPUT GENERIC normalized.tsv
# OUTPUT normalized.tsv
# PARAMETER genome.build [NCBI36, NCBI35] DEFAULT NCBI36 (The genome build to use for fetching the CNV data.)

# count-overlapping-cnvs.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-03-11

dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('This script can only be run on files that have the following columns: chromosome, start, end.')

# load cnvs
if (genome.build=='NCBI35') {
  cnv <- read.table('http://projects.tcag.ca/variation/downloads/variation.hg17.v8.aug.2009.txt', header=TRUE, sep='\t', as.is=TRUE)
} else {
  cnv <- read.table('http://projects.tcag.ca/variation/downloads/variation.hg18.v8.aug.2009.txt', header=TRUE, sep='\t', as.is=TRUE)
}
cnv <- cnv[cnv$VariationType == 'CopyNumber',]
cnv$Chr <- substr(cnv$Chr, 4, 5)

dat$cnvs <- 0
for (i in rownames(dat)) {
  dat[i, 'cnvs'] <- nrow(subset(cnv, Chr == dat[i, 'chromosome'] & 
                                   Start <= dat[i, 'end'] & 
                                     End >= dat[i, 'start']))
}

write.table(dat, file='normalized.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF