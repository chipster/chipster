# TOOL acgh-normalize.R: "Normalize copy number data" (Normalize copy number data before segmentation.)
# INPUT raw.tsv: raw.tsv TYPE GENE_EXPRS 
# OUTPUT normalized.tsv: normalized.tsv
# PARAMETER normalization: normalization TYPE [median: median, mode: mode, none: none] DEFAULT median (Normalization method.)
# PARAMETER number.of.chromosomes: number.of.chromosomes TYPE INTEGER DEFAULT 23 (Number of chromosomes. Usually 23 for sex-matched reference samples and 22 otherwise.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-12-19

source(file.path(chipster.common.path, 'CGHcallPlus.R'))

file <- 'raw.tsv'
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, as.is=TRUE, check.names=FALSE)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

anno <- dat[,-grep("^chip\\.", colnames(dat))]
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
dat2$chromosome <- as.integer(dat2$chromosome)

dat2 <- dat2[order(dat2$chromosome, dat2$start),]

cgh.raw <- make_cghRaw(dat2)
cgh.pre <- preprocess(cgh.raw, nchrom=number.of.chromosomes)
cgh.nor <- normalize(cgh.pre, method=normalization)

dat3 <- cbind(anno[featureNames(cgh.nor),], copynumber(cgh.nor))
colnames(dat3) <- colnames(dat)

dat3$chromosome <- as.character(dat3$chromosome)
dat3$chromosome[dat3$chromosome=='23'] <- 'X'
dat3$chromosome[dat3$chromosome=='24'] <- 'Y'

options(scipen=10)
write.table(dat3, file='normalized.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF
