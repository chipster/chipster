# ANALYSIS "aCGH tools (beta testing)"/"Fetch probe positions from CanGEM" (Fetches microarray probe positions from the CanGEM database.)
# INPUT GENE_EXPRS normalized.tsv
# OUTPUT normalized.tsv
# PARAMETER platform.accession STRING DEFAULT CG-PLM- (The accession of the platform.)
# PARAMETER genome.build [GRCh37, NCBI36, NCBI35] DEFAULT GRCh37 (The genome build to use for adding the chromosome names and start and end base pair positions for the probes.)

# fetch-probe-positions-from-cangem.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-03-18

# check for valid platform accession
platform.accession <- toupper(platform.accession)
if (length(grep('^CG-PLM-[0-9]+$', platform.accession))==0)
	stop('Not a valid platform accession: ', platform.accession)

dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

# remove probe positions if already present
dat$chromosome <- NULL
dat$start <- NULL
dat$end <- NULL
dat$cytoband <- NULL

# load platform
platform <- read.table(paste('http://www.cangem.org/download.php?platform=', platform.accession, '&flag', genome.build, sep=''), sep='\t', header=TRUE, as.is=TRUE)
colnames(platform) <- tolower(colnames(platform))
colnames(platform)[colnames(platform)=='chr'] <- 'chromosome'
rownames(platform) <- platform[,1]
platform$chromosome <- factor(platform$chromosome, levels=c(1:22, "X", "Y", "MT"), ordered=TRUE)

dat2 <- cbind(platform[rownames(dat), c('chromosome', 'start', 'end')], dat, row.names=rownames(dat))
dat2 <- dat2[order(dat2$chromosome, dat2$start),]

write.table(dat2, file='normalized.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF
