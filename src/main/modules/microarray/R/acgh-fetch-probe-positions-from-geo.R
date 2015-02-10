# TOOL acgh-fetch-probe-positions-from-geo.R: "Fetch probe positions from GEO" (Fetches microarray probe positions from the GEO database.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS
# OUTPUT probe-positions.tsv: probe-positions.tsv
# PARAMETER platform: Platform TYPE STRING DEFAULT GPL (The accession of the platform.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-22

source(file.path(chipster.common.path, 'library-Chipster.R'))

# check for valid accession
platform <- toupper(platform)
if (length(grep('^GPL[0-9]+$', platform)) == 0)
  stop('CHIPSTER-NOTE: Not a valid accession: ', platform)

dat <- readData("normalized.tsv")

# remove probe positions if already present
dat$chromosome <- NULL
dat$start <- NULL
dat$end <- NULL
dat$cytoband <- NULL

# load data
library(GEOquery)

gds <- getGEO(platform)
plat <- gds@dataTable@table
plat <- plat[!is.na(plat[, 1]) & plat[, 1] != "", ]
rownames(plat) <- plat[, 1]
colnames(plat) <- toupper(colnames(plat))

# Agilent
if ('CHROMOSOMAL_LOCATION' %in% colnames(plat)) {
  plat$chromosome <- gsub('chr|_random|_hla_hap1|_hla_hap2|:.*','', plat$CHROMOSOMAL_LOCATION)
  plat$start <- as.integer(gsub('.*:|-.*','', plat$CHROMOSOMAL_LOCATION))
  plat$end <- as.integer(gsub('.*-|','', plat$CHROMOSOMAL_LOCATION))
}
if (all(is.na(plat$chromosome)) && 'SYSTEMATICNAME' %in% colnames(plat)) {
  plat$SYSTEMATICNAME[grep('^.*:[0-9]*-[0-9]*$', plat$SYSTEMATICNAME, invert=TRUE)] <- ''
  plat$chromosome <- gsub('chr|_random|_hla_hap1|_hla_hap2|:.*','', plat$SYSTEMATICNAME)
  plat$start <- as.integer(gsub('.*:|-.*','', plat$SYSTEMATICNAME))
  plat$end <- as.integer(gsub('.*-|','', plat$SYSTEMATICNAME))
}
if (all(is.na(plat$cytoband)) && 'CYTOBAND' %in% colnames(plat))
  plat$cytoband <- gsub('hs\\|', '', plat$CYTOBAND)
if (all(is.na(plat$symbol)) && 'GENE_SYMBOL' %in% colnames(plat))
  plat$symbol <- plat$GENE_SYMBOL
if (all(is.na(plat$description)) && 'GENE_NAME' %in% colnames(plat))
  plat$description <- plat$GENE_NAME
if (all(is.na(plat$description)) && 'DESCRIPTION' %in% colnames(plat))
  plat$description <- plat$DESCRIPTION

# Nimblegen
if (all(is.na(plat$chromosome)) && 'CHROMOSOME' %in% colnames(plat))
  plat$chromosome <- plat$CHROMOSOME
if (all(is.na(plat$start)) && 'RANGE_START' %in% colnames(plat))
  plat$start <- as.integer(plat$RANGE_START)
if (all(is.na(plat$end)) && 'RANGE_END' %in% colnames(plat))
  plat$end <- as.integer(plat$RANGE_END)

# other
if (all(is.na(plat$chromosome)) && 'CHROMOSOME_NR' %in% colnames(plat))
  plat$chromosome <- plat$CHROMOSOME_NR
if (all(is.na(plat$start)) && 'START' %in% colnames(plat))
  plat$start <- as.integer(plat$START)
if (all(is.na(plat$start)) && 'POSITION' %in% colnames(plat))
  plat$start <- as.integer(plat$POSITION)
if (all(is.na(plat$start)) && 'KB POSITION' %in% colnames(plat))
  plat$start <- as.integer(plat$'KB POSITION') * 1000
if (all(is.na(plat$end)) && 'END' %in% colnames(plat))
  plat$end <- as.integer(plat$END)
if (all(is.na(plat$symbol)) && 'SYMBOL' %in% colnames(plat))
  plat$symbol <- plat$SYMBOL
if (all(is.na(plat$description)) && 'GENE_DESCRIPTION' %in% colnames(plat))
  plat$description <- plat$GENE_DESCRIPTION

# if missing, impute end column from start (assuming 60 bp probes)
if (all(is.na(plat$end)))
  plat$end <- plat$start + 60

dat2 <- cbind(plat[rownames(dat), c('chromosome', 'start', 'end', 'cytoband', 'symbol', 'description')], dat, row.names=rownames(dat))

# remove empty annotation columns and clean up
for (x in c('chromosome', 'start', 'end', 'cytoband', 'symbol', 'description')) {
  if (all(is.na(dat2[,x]))) {
    dat2[, x] <- NULL
  } else {
    dat2[, x] <- gsub('"|\'|#|\t', '', dat2[, x])
  }
}

writeData(dat2, "probe-positions.tsv")

# EOF
