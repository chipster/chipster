# TOOL acgh-fetch-probe-positions-from-geo.R: "Fetch probe positions from GEO" (Fetches microarray probe positions from the GEO database.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS
# OUTPUT probe-positions.tsv: probe-positions.tsv
# PARAMETER platform: platform TYPE STRING DEFAULT GPL (The accession of the platform.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-12-11

# check for valid accession
platform <- toupper(platform)
if (length(grep('^GPL[0-9]+$', platform)) == 0)
  stop('CHIPSTER-NOTE: Not a valid accession: ', platform)

file <- 'normalized.tsv'
dat2 <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, as.is=TRUE, check.names=FALSE)

# remove probe positions if already present
dat2$chromosome <- NULL
dat2$start <- NULL
dat2$end <- NULL
dat2$cytoband <- NULL

# load data
library(GEOquery)

gds <- getGEO(platform)

dat <- gds@dataTable@table
dat$chromosome <- NA
dat$start <- NA
dat$end <- NA
dat$cytoband <- NA
dat$symbol <- NA
dat$description <- NA

# annotation columns
plat <- gds@dataTable@table
colnames(plat) <- toupper(colnames(plat))

# Agilent
if ('CHROMOSOMAL_LOCATION' %in% colnames(plat)) {
  dat$chromosome <- gsub('chr|_random|_hla_hap1|_hla_hap2|:.*','', plat$CHROMOSOMAL_LOCATION)
  dat$start <- as.integer(gsub('.*:|-.*','', plat$CHROMOSOMAL_LOCATION))
  dat$end <- as.integer(gsub('.*-|','', plat$CHROMOSOMAL_LOCATION))
}
if (all(is.na(dat$cytoband)) && 'CYTOBAND' %in% colnames(plat))
  dat$cytoband <- gsub('hs\\|', '', plat$CYTOBAND)
if (all(is.na(dat$symbol)) && 'GENE_SYMBOL' %in% colnames(plat))
  dat$symbol <- plat$GENE_SYMBOL
if (all(is.na(dat$description)) && 'GENE_NAME' %in% colnames(plat))
  dat$description <- plat$GENE_NAME
if (all(is.na(dat$description)) && 'DESCRIPTION' %in% colnames(plat))
  dat$description <- plat$DESCRIPTION

# Nimblegen
if (all(is.na(dat$chromosome)) && 'CHROMOSOME' %in% colnames(plat))
  dat$chromosome <- plat$CHROMOSOME
if (all(is.na(dat$start)) && 'RANGE_START' %in% colnames(plat))
  dat$start <- as.integer(plat$RANGE_START)
if (all(is.na(dat$end)) && 'RANGE_END' %in% colnames(plat))
  dat$end <- as.integer(plat$RANGE_END)

# other
if (all(is.na(dat$chromosome)) && 'CHROMOSOME_NR' %in% colnames(plat))
  dat$chromosome <- plat$CHROMOSOME_NR
if (all(is.na(dat$start)) && 'START' %in% colnames(plat))
  dat$start <- as.integer(plat$START)
if (all(is.na(dat$end)) && 'END' %in% colnames(plat))
  dat$end <- as.integer(plat$END)
if (all(is.na(dat$symbol)) && 'SYMBOL' %in% colnames(plat))
  dat$symbol <- plat$SYMBOL
if (all(is.na(dat$description)) && 'GENE_DESCRIPTION' %in% colnames(plat))
  dat$description <- plat$GENE_DESCRIPTION

dat3 <- cbind(dat[rownames(dat2), c('chromosome', 'start', 'end', 'cytoband', 'symbol', 'description')], dat2, row.names=rownames(dat2))

# remove empty annotation columns and clean up
for (x in c('chromosome', 'start', 'end', 'cytoband', 'symbol', 'description')) {
  if (all(is.na(dat3[,x]))) {
    dat3[,x] <- NULL
  } else {
    dat3[,x] <- gsub('"|\'|#|\t', '', dat3[,x])
  }
}

options(scipen=10)
write.table(dat3, file='probe-positions.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF
