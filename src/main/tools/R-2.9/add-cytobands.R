# ANALYSIS "aCGH tools (beta testing)"/"Add cytogenetic bands" (Adds the cytogenetic band information using chromosome names and start/end base pair positions present in the data.)
# INPUT GENERIC normalized.tsv
# OUTPUT cytobands.tsv
# PARAMETER genome.build [GRCh37, NCBI36, NCBI35, NCBI34] DEFAULT GRCh37 (The genome build to use. GRCh37 = hg19, NCBI36 = hg18, NCBI35 = hg17, NCBI34 = hg16.)

# add-cytobands.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-09-04

dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

# load cytobands
bands <- read.table(paste('http://www.cangem.org/download.php?platform=CG-PLM-6&flag=', genome.build, sep=''), sep='\t', header=TRUE, as.is=TRUE)
colnames(bands) <- tolower(colnames(bands))
colnames(bands)[colnames(bands)=='chr'] <- 'chromosome'
rownames(bands) <- bands[,1]

dat2 <- dat[,pos]
dat2$cytoband <- NA
dat2 <- cbind(dat2, dat[,setdiff(colnames(dat), pos)])
dat <- dat2

dat.na <- dat[is.na(dat$chromosome),]
dat <- dat[!is.na(dat$chromosome),]

for (band in rownames(bands)) {
  index <- dat$chromosome == bands[band, 'chromosome'] &
                dat$start >= bands[band, 'start'] &
                dat$start <= bands[band, 'end']
  if (length(index)>0)
    dat[index, 'startband'] <- bands[band, 'band']
  index <- dat$chromosome == bands[band, 'chromosome'] &
                  dat$end >= bands[band, 'start'] &
                  dat$end <= bands[band, 'end']
  if (length(index)>0)
    dat[index, 'endband'] <- bands[band, 'band']
}

dat[!is.na(dat$startband), 'cytoband'] <- paste(dat[!is.na(dat$startband), 'startband'], '-', dat[!is.na(dat$startband), 'endband'], sep='')
dat[!is.na(dat$startband) & dat$startband==dat$endband, 'cytoband'] <- dat[!is.na(dat$startband) & dat$startband==dat$endband, 'startband']

dat$startband <- NULL
dat$endband <- NULL

dat <- rbind(dat, dat.na)

write.table(dat, file='cytobands.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF