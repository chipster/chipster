# TOOL cna-segment.R: "Segment copy number data" (Segment copy number data before calling aberrations.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT segmented.tsv: segmented.tsv 
# PARAMETER normalization: normalization TYPE [median: median, mode: mode, none: none] DEFAULT none (Normalization method.)
# PARAMETER number.of.chromosomes: "number of chromosomes" TYPE INTEGER DEFAULT 23 (Number of chromosomes. Usually 23 for sex-matched reference samples and 22 otherwise.)
# PARAMETER minimum.number.of.probes.per.segment: "minimum number of bins per segment" TYPE [2: 2, 3: 3, 4: 4, 5: 5] DEFAULT 2 (Minimum number of bins per segment.)
# PARAMETER minimum.number.of.sds.between.segments: "minimum number of standard deviations between segments" TYPE DECIMAL FROM 0 TO 10 DEFAULT 0 (Minimum number of standard deviations required between segments.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2013-02-23

source(file.path(chipster.common.path, 'CGHcallPlus.R'))

file <- 'normalized.tsv'
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, as.is=TRUE, check.names=FALSE)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

logratios <- dat[,grep("^chip\\.", colnames(dat))]
logratios[logratios == -Inf] <- -10

dat2 <- data.frame(probe=rownames(dat), dat[,c('chromosome', 'start', 'end')], logratios, stringsAsFactors=FALSE)
chips <- grep("^chip\\.", colnames(dat), value=TRUE)

if (ncol(dat2)==4)
  stop('CHIPSTER-NOTE: No array data found. The input file must have columns labeled with "chip.".')
if (ncol(dat2)==5)
  colnames(dat2)[5] <- chips[1]

dat2 <- dat2[!is.na(dat2$chromosome),]
dat2 <- dat2[!is.na(dat2$start),]
dat2 <- dat2[!is.na(dat2$end),]

dat2$chromosome[dat2$chromosome=='X'] <- 23
dat2$chromosome[dat2$chromosome=='Y'] <- 24
dat2$chromosome[dat2$chromosome=='MT'] <- 25
dat2$chromosome <- as.integer(dat2$chromosome)

dat2 <- dat2[order(dat2$chromosome, dat2$start),]

cgh.raw <- make_cghRaw(dat2)
cgh.pre <- preprocess(cgh.raw, nchrom=number.of.chromosomes)
cgh.nor <- normalize(cgh.pre, method=normalization)
cgh.seg <- segmentData(cgh.nor, min.width=as.integer(minimum.number.of.probes.per.segment), undo.splits='sdundo', undo.SD=minimum.number.of.sds.between.segments)

dat3 <- data.frame(cgh.seg@featureData@data)
colnames(dat3) <- c('chromosome', 'start', 'end')

for (col in c('cytoband', 'symbol', 'description', 'cnvs'))
  if (col %in% colnames(dat))
    dat3[,col] <- dat[rownames(dat3), col]

copynumber <- copynumber(cgh.seg)
colnames(copynumber) <- chips

dat3 <- cbind(dat3, copynumber)

segmented <- segmented(cgh.seg)
colnames(segmented) <- sub('chip.', 'segmented.', chips)
dat3 <- cbind(dat3, segmented)

dat3$chromosome <- as.character(dat3$chromosome)
dat3$chromosome[dat3$chromosome=='23'] <- 'X'
dat3$chromosome[dat3$chromosome=='24'] <- 'Y'
dat3$chromosome[dat3$chromosome=='25'] <- 'MT'

options(scipen=10)
write.table(dat3, file='segmented.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF
