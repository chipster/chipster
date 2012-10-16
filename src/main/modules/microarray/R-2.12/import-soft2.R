# TOOL import-soft2.R: "Import from GEO" (Import a data set directly from GEO. Be sure to specify the chiptype as an Affymetrix chip name, or either Illumina or cDNA.)
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER GDS.name: accession TYPE STRING DEFAULT GDS858 (GDS or GSE number of the experiment.)
# PARAMETER platform: platform TYPE STRING DEFAULT GPL (In case the series contains multiple platforms, specify the accession of the platform to import. If there is just one, this platform is ignored.)
# PARAMETER chiptype: chiptype TYPE STRING DEFAULT cDNA (If the microarray platform used is an Affymetrix one, the name of the Bioconductor annotation package, Illumina for Illumina arrays, or cDNA for everything else.)

# 2012-09-26
# Ilari Scheinin <firstname.lastname@gmail.com>

# JTT 9.8.2007

# check for valid accession
GDS.name <- toupper(GDS.name)
if (length(grep('^(GDS|GSE|)[0-9]+$', GDS.name)) == 0)
  stop('CHIPSTER-NOTE: Not a valid accession: ', GDS.name)

# load data
library(GEOquery)

gds <- getGEO(GDS.name)
if (class(gds) == 'GDS') {
  eset <- GDS2eSet(gds)
} else if (length(gds) == 1) {
  eset <- gds[[1]]
} else {
  w <- grep(platform, names(gds))
  if (length(w) != 1)
    stop('CHIPSTER-NOTE: Please use the platform argument (currently set to "', platform, '") to specify which data set to load. Available data sets:\n', paste(names(gds), collapse='\n')) 
  eset <- gds[[w]]
}

# generate phenodata
if ('geo_accession' %in% colnames(pData(eset))) {
  sample <- pData(eset)$geo_accession
} else if ('sample' %in% colnames(pData(eset))) {
  sample <- pData(eset)$sample
} else {
  sample <- sprintf('microarray%.3i', 1:ncol(eset))
}
group <- rep('', length(sample))
phenodata <- data.frame(sample=sample, original_name=sample, chiptype=chiptype, group=group, description=sample, pData(eset))

dat <- data.frame(chromosome=NA, start=NA, end=NA, cytoband=NA, symbol=NA, description=NA, exprs(eset))
colnames(dat)[-(1:6)] <- paste('chip.', sample, sep='')

# GDS annotation columns
if ('Gene symbol' %in% colnames(eset@featureData@data))
  dat$symbol <- eset@featureData@data[, 'Gene symbol']
if ('Gene title' %in% colnames(eset@featureData@data))
  dat$description <- eset@featureData@data[, 'Gene title']

# Agilent annotation columns
if ('CHROMOSOMAL_LOCATION' %in% colnames(eset@featureData@data)) {
  dat$chromosome <- gsub('chr|_random|_hla_hap1|_hla_hap2|:.*','', eset@featureData@data$CHROMOSOMAL_LOCATION)
  dat$start <- as.integer(gsub('.*:|-.*','', eset@featureData@data$CHROMOSOMAL_LOCATION))
  dat$end <- as.integer(gsub('.*-|','', eset@featureData@data$CHROMOSOMAL_LOCATION))
}
if (all(is.na(dat$cytoband)) && 'CYTOBAND' %in% colnames(eset@featureData@data))
  dat$cytoband <- gsub('hs\\|', '', eset@featureData@data$CYTOBAND)
if (all(is.na(dat$symbol)) && 'GENE_SYMBOL' %in% colnames(eset@featureData@data))
  dat$symbol <- eset@featureData@data$GENE_SYMBOL
if (all(is.na(dat$description)) && 'GENE_NAME' %in% colnames(eset@featureData@data))
  dat$description <- eset@featureData@data$GENE_NAME
if (all(is.na(dat$description)) && 'DESCRIPTION' %in% colnames(eset@featureData@data))
  dat$description <- eset@featureData@data$DESCRIPTION

# Nimblegen annotation columns
if (all(is.na(dat$chromosome)) && 'CHROMOSOME' %in% colnames(eset@featureData@data))
  dat$chromosome <- eset@featureData@data$CHROMOSOME
if (all(is.na(dat$start)) && 'RANGE_START' %in% colnames(eset@featureData@data))
  dat$start <- as.integer(eset@featureData@data$RANGE_START)
if (all(is.na(dat$end)) && 'RANGE_END' %in% colnames(eset@featureData@data))
  dat$end <- as.integer(eset@featureData@data$RANGE_END)

# remove empty annotation columns
for (x in c('chromosome', 'start', 'end', 'cytoband', 'symbol', 'description'))
  if (all(is.na(dat[,x])))
    dat[,x] <- NULL

# write output files
options(scipen=10)
write.table(dat, file='normalized.tsv', quote=FALSE, sep='\t')
write.table(phenodata, file='phenodata.tsv', quote=FALSE, sep='\t', row.names=FALSE)

# EOF
