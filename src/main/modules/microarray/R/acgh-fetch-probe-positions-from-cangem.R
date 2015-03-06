# TOOL acgh-fetch-probe-positions-from-cangem.R: "Fetch probe positions from CanGEM" (Fetches microarray probe positions from the CanGEM database.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT probe-positions.tsv: probe-positions.tsv 
# PARAMETER platform.accession: "Platform accession" TYPE STRING DEFAULT CG-PLM- (The accession of the platform.)
# PARAMETER genome.build: "Genome build" TYPE [GRCh37: GRCh37, NCBI36: NCBI36, NCBI35: NCBI35, NCBI34: NCBI34] DEFAULT GRCh37 (The genome build to use for adding the chromosome names and start and end base pair positions for the probes.)
# PARAMETER username: Username TYPE STRING DEFAULT empty (Username, in case the data is password-protected. WARNING: This will store your username password in the Chipster history files. To avoid this, use the session parameter.)
# PARAMETER password: Password TYPE STRING DEFAULT empty (Password, in case the data is password-protected. WARNING: This will store your username password in the Chipster history files. To avoid this, use the session parameter.)
# PARAMETER session: "Session ID" TYPE STRING DEFAULT empty (Session ID. To avoid saving your username password in Chipster history files, log in at http: www.cangem.org using a web browser, then copy&paste your session ID from the lower right corner of the CanGEM website. This will allow Chipster to access your password-protected data until you log out of the web site (or the session times out\).)

# fetch-probe-positions-from-cangem.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-28

source(file.path(chipster.common.path, 'library-Chipster.R'))

# check for valid platform accession
platform.accession <- toupper(platform.accession)
if (length(grep('^CG-PLM-[0-9]+$', platform.accession)) == 0)
  stop('CHIPSTER-NOTE: Not a valid platform accession: ', platform.accession)

dat <- readData("normalized.tsv")

# remove probe positions if already present
dat$chromosome <- NULL
dat$start <- NULL
dat$end <- NULL
dat$cytoband <- NULL

# construct the string used in authenticating
if (session != 'empty' && session != '') {
  auth <- paste('&PHPSESSID=', session, sep='')
} else if (username != 'empty' && username != '' && password != 'empty' && password != '') {
  auth <- paste('&username=', username, '&password=', password, sep = '')
} else auth <- ''

# load platform
plat <- read.table(paste('http://www.cangem.org/download.php?platform=', platform.accession, '&flag=', genome.build, auth, sep=''), sep='\t', header=TRUE, as.is=TRUE)
colnames(plat) <- tolower(colnames(plat))
colnames(plat)[colnames(plat) == 'chr'] <- 'chromosome'
rownames(plat) <- plat[, 1]

dat2 <- cbind(plat[rownames(dat), c('chromosome', 'start', 'end')], dat, row.names=rownames(dat))
dat2 <- dat2[order(chromosomeToInteger(dat2$chromosome), dat2$start), ]

writeData(dat2, "probe-positions.tsv")

# EOF
