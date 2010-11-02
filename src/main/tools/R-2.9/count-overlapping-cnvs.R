# ANALYSIS "aCGH tools"/"Count overlapping CNVs" (Counts overlapping CNVs from the database of genomic variants.)
# INPUT GENERIC normalized.tsv
# OUTPUT cnvs.tsv
# PARAMETER genome.build [GRCh37, NCBI36, NCBI35] DEFAULT GRCh37 (The genome build to use. GRCh37 = hg19, NCBI36 = hg18, NCBI35 = hg17, NCBI34 = hg16.)

# count-overlapping-cnvs.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2010-10-05

dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

dat$chromosome <- factor(dat$chromosome, levels=c(1:22, 'X', 'Y', 'MT'), ordered=TRUE)

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
cnv$Chr <- factor(substr(cnv$Chr, 4, 5), levels=c(1:22, 'X', 'Y', 'MT'), ordered=TRUE)

cnv <- cnv[order(cnv$Chr, cnv$Start, cnv$End), c('Chr', 'Start', 'End')]
joined <- data.frame()
prev <- cnv[1,]
for (i in 2:nrow(cnv)) {
  if (cnv[i, 'Chr'] != prev$Chr || cnv[i, 'Start'] > (prev$End + 1)) {
    joined <- rbind(joined, prev)
    prev <- cnv[i,]
  } else
    prev$End <- max(prev$End, cnv[i, 'End'])
}
joined <- rbind(joined, prev)

for (i in rownames(dat2)) {
  chr <- dat2[i, 'chromosome']
  start <- dat2[i, 'start']
  end <- dat2[i, 'end']
  dat2[i, 'cnv.count'] <- nrow(cnv[cnv$Chr   == chr &
                                   cnv$Start <= end &
                                   cnv$End   >= start,])
  overlaps <- joined[joined$Chr   == chr &
                     joined$Start <= end &
                     joined$End   >= start,]
  bases <- 0
  for (j in rownames(overlaps))
    bases <- bases + min(end, overlaps[j, 'End']) - max(start, overlaps[j, 'Start']) + 1
  dat2[i, 'cnv.per.Mb'] <- round(bases / (end - start + 1) * 1000000)
}

dat2 <- cbind(dat2, dat[,first.data.col:ncol(dat)])

write.table(dat2, file='cnvs.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF