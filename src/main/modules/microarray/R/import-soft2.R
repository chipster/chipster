# TOOL import-soft2.R: "Import from GEO" (Import a data set directly from GEO. Be sure to specify the chiptype as an Affymetrix chip name, or either Illumina or cDNA.)
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER GDS.name: Accession TYPE STRING DEFAULT GSE (GDS or GSE number of the experiment.)
# PARAMETER platform: Platform TYPE STRING DEFAULT GPL (In case the series contains multiple platforms, specify the accession of the platform to import. If there is just one, this platform is ignored.)
# PARAMETER chiptype: "Affymetrix/Illumina chiptype" TYPE STRING DEFAULT other (If the microarray platform used is an Affymetrix one, the name of the Bioconductor annotation package. For Illumina arrays, fill in Illumina. For everything else, either cDNA or other.)
# PARAMETER matrix: "Matrix number" TYPE INTEGER FROM 1 TO 100 DEFAULT 1 (The rank number of softmatrix to retrieve. Please note that most GEO-entries include only one softmatrix)

# JTT 09.08.2007
# IS: 03.18.2013: Ilari Scheinin <firstname.lastname@gmail.com>
# MK: 11.12.2013: Added selector for specifying which dataset to load, if the entry has multiple one with the same GPL identifier

# check for valid accession
GDS.name <- toupper(GDS.name)
if (length(grep('^(GDS|GSE|)[0-9]+$', GDS.name)) == 0)
  stop('CHIPSTER-NOTE: Not a valid accession: ', GDS.name)

if (chiptype=='other')
  chiptype <- 'cDNA'

# load data
library(GEOquery)

gds <- getGEO(GDS.name)
if (class(gds) == 'GDS') {
  eset <- GDS2eSet(gds)
} else if (length(gds) == 1) {
  eset <- gds[[1]]
} else {
  w <- grep(platform, names(gds))
  if (length(w) != 1) {
    eset <- gds[[w[matrix]]]
  } else {
    eset <- gds[[w]]
  }
}

# clean up phenodata
pdata <- pData(eset)
for (x in colnames(pdata))
  pdata[,x] <- gsub('"|\'|#|\t', '', pdata[,x])

# generate phenodata
if ('geo_accession' %in% colnames(pdata)) {
  sample <- pdata$geo_accession
} else if ('sample' %in% colnames(pdata)) {
  sample <- pdata$sample
} else {
  sample <- sprintf('microarray%.3i', 1:ncol(eset))
}
group <- rep('', length(sample))
phenodata <- data.frame(sample=sample, original_name=sample, chiptype=chiptype, group=group, description=sample, pdata)

dat <- data.frame(chromosome=NA, start=NA, end=NA, cytoband=NA, symbol=NA, description=NA, exprs(eset))
colnames(dat)[-(1:6)] <- paste('chip.', sample, sep='')

# annotation columns
plat <- eset@featureData@data
colnames(plat) <- toupper(colnames(plat))

# GDS
if ('Gene symbol' %in% colnames(plat))
  dat$symbol <- plat[, 'Gene symbol']
if ('Gene title' %in% colnames(plat))
  dat$description <- plat[, 'Gene title']

# Agilent
if ('CHROMOSOMAL_LOCATION' %in% colnames(plat)) {
  dat$chromosome <- gsub('chr|_random|_hla_hap1|_hla_hap2|:.*','', plat$CHROMOSOMAL_LOCATION)
  dat$start <- as.integer(gsub('.*:|-.*','', plat$CHROMOSOMAL_LOCATION))
  dat$end <- as.integer(gsub('.*-|','', plat$CHROMOSOMAL_LOCATION))
}
if (all(is.na(dat$chromosome)) && 'SYSTEMATICNAME' %in% colnames(plat)) {
  plat$SYSTEMATICNAME[grep('^.*:[0-9]*-[0-9]*$', plat$SYSTEMATICNAME, invert=TRUE)] <- ''
  dat$chromosome <- gsub('chr|_random|_hla_hap1|_hla_hap2|:.*','', plat$SYSTEMATICNAME)
  dat$start <- as.integer(gsub('.*:|-.*','', plat$SYSTEMATICNAME))
  dat$end <- as.integer(gsub('.*-|','', plat$SYSTEMATICNAME))
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
if (all(is.na(dat$start)) && 'POSITION' %in% colnames(plat))
  dat$start <- as.integer(plat$POSITION)
if (all(is.na(dat$start)) && 'KB POSITION' %in% colnames(plat))
  dat$start <- as.integer(plat$'KB POSITION') * 1000
if (all(is.na(dat$end)) && 'END' %in% colnames(plat))
  dat$end <- as.integer(plat$END)
if (all(is.na(dat$symbol)) && 'SYMBOL' %in% colnames(plat))
  dat$symbol <- plat$SYMBOL
if (all(is.na(dat$description)) && 'GENE_DESCRIPTION' %in% colnames(plat))
  dat$description <- plat$GENE_DESCRIPTION

# if missing, impute end column from start (assuming 60 bp probes)
if (all(is.na(dat$end)))
  dat$end <- dat$start + 60

# remove empty annotation columns and clean up
for (x in c('chromosome', 'start', 'end', 'cytoband', 'symbol', 'description')) {
  if (all(is.na(dat[,x]))) {
    dat[,x] <- NULL
  } else {
    dat[,x] <- gsub('"|\'|#|\t', '', dat[,x])
  }
}

# write output files
options(scipen=10)
write.table(dat, file='normalized.tsv', quote=FALSE, sep='\t')
write.table(phenodata, file='phenodata.tsv', quote=FALSE, sep='\t', row.names=FALSE)

# EOF
