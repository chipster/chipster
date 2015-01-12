# TOOL acgh-count-overlapping-cnvs.R: "Count overlapping CNVs" (Counts overlapping CNVs from the database of genomic variants.)
# INPUT normalized.tsv: normalized.tsv TYPE GENERIC 
# OUTPUT cnvs.tsv: cnvs.tsv 
# PARAMETER genome.build: "Genome build" TYPE [GRCh37: GRCh37, NCBI36: NCBI36] DEFAULT GRCh37 (The genome build to use. GRCh37 = hg19, NCBI36 = hg18.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-28

source(file.path(chipster.common.path, 'library-Chipster.R'))

dat <- readData("normalized.tsv")

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

dat$chromosome <- factor(dat$chromosome, levels=c(1:22, 'X', 'Y', 'MT'), ordered=TRUE)

first.data.col <- max(0, min(grep('chip', names(dat)), grep('flag', names(dat))))

if (first.data.col > 0) {
  dat2 <- dat[1:first.data.col-1]
} else {
  dat2 <- dat
}

# load cnvs
if (genome.build == 'NCBI36') {
  cnv <- read.table('http://dgv.tcag.ca/dgv/docs/NCBI36_hg18_variants_2013-07-23.txt', header=TRUE, sep='\t', as.is=TRUE)
} else {
  cnv <- read.table('http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2013-07-23.txt', header=TRUE, sep='\t', as.is=TRUE)
}
cnv <- cnv[cnv$varianttype == 'CNV',]
cnv$chr <- factor(cnv$chr, levels=c(1:22, 'X', 'Y', 'MT'), ordered=TRUE)
cnv <- cnv[!is.na(cnv$chr),]

cnv <- cnv[order(cnv$chr, cnv$start, cnv$end), c('chr', 'start', 'end')]
joined <- data.frame()
prev <- cnv[1,]
for (i in 2:nrow(cnv)) {
  if (cnv[i, 'chr'] != prev$chr || cnv[i, 'start'] > (prev$end + 1)) {
    joined <- rbind(joined, prev)
    prev <- cnv[i,]
  } else
    prev$end <- max(prev$end, cnv[i, 'end'])
}
joined <- rbind(joined, prev)

cnv.counter <- function(x) {
  chr <- x['chromosome']
  start <- as.integer(x['start'])
  end <- as.integer(x['end'])

  count <- nrow(cnv[cnv$chr   == chr &
                    cnv$start <= end &
                    cnv$end   >= start,])
  overlaps <- joined[joined$chr   == chr &
                     joined$start <= end &
                     joined$end   >= start,]
  bases <- 0
  for (j in rownames(overlaps))
    bases <- bases + min(end, overlaps[j, 'end']) - max(start, overlaps[j, 'start']) + 1
  c(count, signif(bases / (end - start + 1), digits=3))
}

# first try parallel computing
prob <- TRUE
try({
  library(snowfall)
  sfInit(parallel=TRUE, cpus=4)
  sfExport(list=c('cnv', 'joined'))
  dat2[,c('cnv.count', 'cnv.proportion')] <- t(sfApply(dat2, 1, cnv.counter))
  sfStop()
  prob <- FALSE
}, silent=TRUE)
# if problems, fall back to sequential computing
if (prob)
  dat2[,c('cnv.count', 'cnv.proportion')] <- t(apply(dat2, 1, cnv.counter))

if (first.data.col > 0)
  dat2 <- cbind(dat2, dat[,first.data.col:ncol(dat)])

writeData(dat2, "cnvs.tsv")

# EOF
