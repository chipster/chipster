# TOOL acgh-add-cytobands.R: "Add cytogenetic bands" (Adds the cytogenetic band information using chromosome names and start end base pair positions. If this position information is not present in your data set, please first run the Fetch probe positions from GEO/CanGEM tool.)
# INPUT normalized.tsv: normalized.tsv TYPE GENERIC 
# OUTPUT cytobands.tsv: cytobands.tsv 
# PARAMETER genome.build: "Genome build" TYPE [GRCh37: GRCh37, NCBI36: NCBI36] DEFAULT GRCh37 (The genome build to use. GRCh37 = hg19, NCBI36 = hg18.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-08-18

source(file.path(chipster.common.path, 'library-Chipster.R'))

dat <- readData("normalized.tsv")

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat))) != 0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

# load cytobands
cytodir <- list.files(file.path(chipster.tools.path, "genomes", "genomebrowser", "Homo_sapiens"), pattern=paste0("^", genome.build, "\\..*$"), full.name=TRUE)
cytofile <- file.path(cytodir, "cytoband-chr.txt")

if (length(cytofile) == 0)
  stop('CHIPSTER-NOTE: Cytoband file not found.')

if (length(cytofile) > 1)
  cytoband <- sort(cytoband, decreasing=TRUE)

bands <- read.table(cytofile[1], sep="\t", as.is=TRUE, col.names=c("chromosome", "index", "start", "end", "band", "dye"))
bands$chromosome <- sub("^chr", "", bands$chromosome)
bands$band <- paste0(bands$chromosome, bands$band)
rownames(bands) <- bands$band
bands <- bands[bands$chromosome %in% dat$chromosome, ]
bands <- bands[order(bands$index), ]

dat2 <- dat[,pos]
dat2$cytoband <- NA
dat2 <- cbind(dat2, dat[,setdiff(colnames(dat), pos)])
dat <- dat2

for (band in rownames(bands)) {
  index <- !is.na(dat$chromosome) &
    dat$chromosome == bands[band, 'chromosome'] &
    dat$start      >= bands[band, 'start'] &
    dat$start      <= bands[band, 'end']
  if (length(index) > 0)
    dat[index, 'startband'] <- bands[band, 'band']
  index <- !is.na(dat$chromosome) &
    dat$chromosome == bands[band, 'chromosome'] &
    dat$end        >= bands[band, 'start'] &
    dat$end        <= bands[band, 'end']
  if (length(index) > 0)
    dat[index, 'endband'] <- bands[band, 'band']
}

dat$startband[is.na(dat$startband)] <- 'unknown'
dat$endband[is.na(dat$endband)] <- 'unknown'

dat$cytoband <- paste(dat$startband, '-', dat$endband, sep='')
dat$cytoband[dat$startband==dat$endband] <- dat$startband[dat$startband==dat$endband]

dat$startband <- NULL
dat$endband <- NULL

writeData(dat, "cytobands.tsv")

# EOF
