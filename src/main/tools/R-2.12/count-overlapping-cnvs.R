# ANALYSIS "aCGH"/"Count overlapping CNVs" (Counts overlapping CNVs from the database of genomic variants.)
# INPUT GENERIC normalized.tsv
# OUTPUT cnvs.tsv
# PARAMETER genome.build [GRCh37, NCBI36, NCBI35] DEFAULT GRCh37 (The genome build to use. GRCh37 = hg19, NCBI36 = hg18, NCBI35 = hg17, NCBI34 = hg16.)

# count-overlapping-cnvs.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-04-06

dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

dat$chromosome <- factor(dat$chromosome, levels=c(1:22, 'X', 'Y', 'MT'), ordered=TRUE)

first.data.col <- min(0, grep('chip', names(dat)), grep('flag', names(dat)))

if (first.data.col > 0) {
  dat2 <- dat[1:first.data.col-1]
} else {
  dat2 <- dat
}

# load cnvs
if (genome.build=='NCBI35') {
  # cnv <- read.table('http://projects.tcag.ca/variation/downloads/variation.hg17.v10.nov.2010.txt', header=TRUE, sep='\t', as.is=TRUE)
  cnv <- read.table(file.path(chipster.tools.path, 'DGV', 'variation.hg17.v10.nov.2010.txt'), header=TRUE, sep='\t', as.is=TRUE)
} else if (genome.build=='NCBI36') {
  # cnv <- read.table('http://projects.tcag.ca/variation/downloads/variation.hg18.v10.nov.2010.txt', header=TRUE, sep='\t', as.is=TRUE)
  cnv <- read.table(file.path(chipster.tools.path, 'DGV', 'variation.hg18.v10.nov.2010.txt'), header=TRUE, sep='\t', as.is=TRUE)
} else {
  # cnv <- read.table('http://projects.tcag.ca/variation/downloads/variation.hg19.v10.nov.2010.txt', header=TRUE, sep='\t', as.is=TRUE)
  cnv <- read.table(file.path(chipster.tools.path, 'DGV', 'variation.hg19.v10.nov.2010.txt'), header=TRUE, sep='\t', as.is=TRUE)
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

cnv.counter <- function(x) {
  chr <- x['chromosome']
  start <- as.integer(x['start'])
  end <- as.integer(x['end'])

  count <- nrow(cnv[cnv$Chr   == chr &
                    cnv$Start <= end &
                    cnv$End   >= start,])
  overlaps <- joined[joined$Chr   == chr &
                     joined$Start <= end &
                     joined$End   >= start,]
  bases <- 0
  for (j in rownames(overlaps))
    bases <- bases + min(end, overlaps[j, 'End']) - max(start, overlaps[j, 'Start']) + 1
  c(count, round(bases / (end - start + 1) * 1000000))
}

# first try parallel computing
prob <- TRUE
try({
  library(snowfall)
  sfInit(parallel=TRUE, cpus=4)
  sfExport(list=c('cnv', 'joined'))
  dat2[,c('cnv.count', 'cnv.per.Mb')] <- t(sfApply(dat2, 1, cnv.counter))
  sfStop()
  prob <- FALSE
}, silent=TRUE)
# if problems, fall back to sequential computing
if (prob)
  dat2[,c('cnv.count', 'cnv.per.Mb')] <- t(apply(dat2, 1, cnv.counter))

if (first.data.col > 0)
  dat2 <- cbind(dat2, dat[,first.data.col:ncol(dat)])

write.table(dat2, file='cnvs.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF
