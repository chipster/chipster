# TOOL import-soft2.R: "Import from GEO" (Import a SOFT-formatted datafile directly from GEO. Be sure to specify the chiptype as an Affymetrix chip name, or either Illumina or cDNA.)
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# PARAMETER GDS.name: GDS.name TYPE STRING DEFAULT GDS858 (GDS or GSE number of the experiment.)
# PARAMETER platform: platform TYPE STRING DEFAULT GPL (In case the series contains multiple platforms, specify the accession of the platform to import. If there is just one, this platform is ignored.)
# PARAMETER chiptype: chiptype TYPE STRING DEFAULT cDNA ()

# 2012-09-26
# Ilari Scheinin <firstname.lastname@gmail.com>

# JTT 9.8.2007

# Loads the libraries
library(GEOquery)

# Loads and parses the SOFT file
if(GDS.name == 'empty' || GDS.name == '')
  stop('CHIPSTER-NOTE: You need to specify a valid GDS/GSE!')

gds <- getGEO(GDS.name)
if (length(gds) == 1) {
  eset <- gds[[1]]
} else {
  w <- grep(platform, names(gds))
  if (length(w) != 1)
    stop('CHIPSTER-NOTE: Please use the platform argument (currently set to "', platform, '") to specify which data set to load. Available data sets:\n', paste(names(gds), collapse='\n')) 
  eset <- gds[[w]]
}

# generate phenodata
sample <- pData(eset)$geo_accession
group <- c(rep('', nrow(pData(eset))))
phenodata <- data.frame(sample=sample, original_name=sample, chiptype=chiptype, group=group, description=sample, pData(eset))

dat <- data.frame(chromosome=NA, start=NA, end=NA, cytoband=NA, symbol=NA, description=NA, exprs(eset))
colnames(dat)[-(1:6)] <- paste('chip.', sample, sep='')

# Agilent annotation columns
if ('CHROMOSOMAL_LOCATION' %in% colnames(eset@featureData@data)) {
  dat$chromosome <- gsub('chr|_random|_hla_hap1|_hla_hap2|:.*','', eset@featureData@data$CHROMOSOMAL_LOCATION)
  dat$start <- as.integer(gsub('.*:|-.*','', eset@featureData@data$CHROMOSOMAL_LOCATION))
  dat$end <- as.integer(gsub('.*-|','', eset@featureData@data$CHROMOSOMAL_LOCATION))
}
if ('CYTOBAND' %in% colnames(eset@featureData@data))
  dat$cytoband <- gsub('hs\\|', '', eset@featureData@data$CYTOBAND)
if ('GENE_SYMBOL' %in% colnames(eset@featureData@data))
  dat$symbol <- eset@featureData@data$GENE_SYMBOL
if ('GENE_NAME' %in% colnames(eset@featureData@data))
  dat$description <- eset@featureData@data$GENE_NAME

# Nimblegen annotation columns
if ('CHROMOSOME' %in% colnames(eset@featureData@data))
  dat$chromosome <- eset@featureData@data$CHROMOSOME
if ('RANGE_START' %in% colnames(eset@featureData@data))
  dat$start <- as.integer(eset@featureData@data$RANGE_START)
if ('RANGE_END' %in% colnames(eset@featureData@data))
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
