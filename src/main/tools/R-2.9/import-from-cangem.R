# ANALYSIS "aCGH tools (beta testing)"/"Import from CanGEM" (Load a microarray data set from the CanGEM database.)
# OUTPUT normalized.tsv, phenodata.tsv
# PARAMETER accession STRING (Accession of either a data set, an experiment, a series, or single microarray results.)
# PARAMETER username STRING DEFAULT empty (Username, in case the data is password-protected. WARNING: This will store your username/password in the Chipster history files. To avoid this, use the session parameter.)
# PARAMETER password STRING DEFAULT empty (Password, in case the data is password-protected. WARNING: This will store your username/password in the Chipster history files. To avoid this, use the session parameter.)
# PARAMETER session STRING DEFAULT empty (Session ID. To avoid saving your username/password in Chipster history files, log in at http://www.cangem.org/ using a web browser, then copy&paste your session ID from the lower right corner of the CanGEM website. This will allow Chipster to access your password-protected data until you log out of the web site (or the session times out).)
# PARAMETER agilent.filtering [yes, no] DEFAULT yes (Whether to filter outliers from Agilent arrays. Will be ignored, if downloaded files are not in Agilent file format. Check the help file for details about the filtering function.)
# PARAMETER intra.array.normalization [none, median, loess] DEFAULT loess (Intra-array normalization method for Agilent arrays. Will be ignored, if downloaded files are not in Agilent file format.)
# PARAMETER inter.array.normalization [none, quantile, scale] DEFAULT none (Inter-array normalization method for Agilent arrays. Will be ignored, if downloaded files are not in Agilent file format.)
# PARAMETER affymetrix.normalization [gcrma, rma, mas5] DEFAULT gcrma (Normalization method for Affymetrix arrays. Will be ignored, if downloaded files are not in Affymetrix file format.)
# PARAMETER genome.build [none, GRCh37, NCBI36, NCBI35, NCBI34] DEFAULT GRCh37 (The genome build to use for adding the chromosome names and start and end base pair positions for the probes.)

# import-from-cangem.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-03-20

# check for valid accession
accession <- toupper(accession)
if (length(grep('^CG-(SET|EXP|SER|RES|SAM|PLM)-[0-9]+$', accession))==0)
	stop('Not a valid accession: ', accession)

# construct the string used in authenticating
if (session != 'empty' && session != '') {
	auth <- paste('&PHPSESSID=', session, sep='')
} else if (username != 'empty' && username != '' && password != 'empty' && password != '') {
	auth <- paste('&username=', username, '&password=', password, sep = '')
} else auth <- ''

# fetch list of arrays from CanGEM
cangem.samples <- read.table(paste('http://www.cangem.org/scripts/listhybs.php?dataset=', accession, auth, sep=''), sep='\t', header=TRUE, as.is=TRUE)

# check that we did get some results
if (nrow(cangem.samples)==0)
	stop('No array data found.')

# check for single platform, array type and file format
if (length(unique(cangem.samples$PlatformAccession)) > 1)
	stop('Multiple platforms used: ', paste(unique(cangem.samples$Platform), collapse=', '))
if (length(unique(cangem.samples$Type)) > 1)
	stop('Multiple array types used: ', paste(unique(cangem.samples$Type), collapse=', '))
if (length(unique(cangem.samples$Format)) > 1)
	stop('Multiple file formats used: ', paste(unique(cangem.samples$Format), collapse=', '))

microarrays <- sprintf('microarray%.3i', 1:nrow(cangem.samples)) # or .cel/.tsv ???
chips <- paste('chip.', microarrays, sep='')
arrays <- list()
onecolor <- FALSE
if (cangem.samples$Format[1] == 'agilent') {
	library(limma)
	# define filtering function
	if (agilent.filtering == 'yes') {
		cangem.filter <- function(x) {
			# required fields and their values, field.name=required.value
			reqs <- c(
					ControlType=0,
					gIsBGNonUnifOL=0,
					gIsBGPopnOL=0,
					gIsFeatNonUnifOL=0,
					gIsFeatPopnOL=0,
					gIsFound=1,
					gIsPosAndSignif=1,
					gIsSaturated=0,
					gIsWellAboveBG=1,
					gSurrogateUsed=0,
					IsManualFlag=0,
					rIsBGNonUnifOL=0,
					rIsBGPopnOL=0,
					rIsFeatNonUnifOL=0,
					rIsFeatPopnOL=0,
					rIsFound=1,
					rIsPosAndSignif=1,
					rIsSaturated=0,
					rIsWellAboveBG=1,
					rSurrogateUsed=0)
			# start with a weight of 1 for every spot
			weights <- rep(1, length=nrow(x))
			# check which fields are present, and check if they have the required value
			for (req in names(reqs))
				if (req %in% colnames(x))
					weights <- weights & x[,req] == reqs[req]
			as.numeric(weights)
		}
	} else {
		cangem.filter <- function(x) { 1 }
	}
	# load arrays one by one to be able to handle with files with different number of rows
	for (i in 1:nrow(cangem.samples)) {
		prob <- TRUE
		if (!onecolor) {
			# first try to read the 'normal' columns (means), but if it fails, try medians
			try({
						array.raw <- read.maimages(paste('http://www.cangem.org/download.php?hybridization=', cangem.samples[i,'Accession'], auth, sep=''),
								source=cangem.samples[i,'Format'], wt.fun=cangem.filter)
						prob <- FALSE
					}, silent=TRUE)
			if (prob) {
				try({
							array.raw <- read.maimages(paste('http://www.cangem.org/download.php?hybridization=', cangem.samples[i,'Accession'], auth, sep=''),
									source=cangem.samples[i,'Format'], columns=list(G='gMedianSignal', Gb='gBGMedianSignal', R='rMedianSignal', Rb='rBGMedianSignal'), wt.fun=cangem.filter)
							prob <- FALSE
						}, silent=TRUE)
			}
		}
		if (prob) {
			try({
						array.raw <- read.maimages(paste('http://www.cangem.org/download.php?hybridization=', cangem.samples[i,'Accession'], auth, sep=''),
								source=cangem.samples[i,'Format'], columns=list(G='gMedianSignal', Gb='gBGMedianSignal', R='gMedianSignal', Rb='gBGMedianSignal'), wt.fun=cangem.filter)
						prob <- FALSE
						onecolor <- TRUE
					}, silent=TRUE)
		}
		if (prob) {
			stop('Could not read file ', cangem.samples[i, 'Name'])
		}
		array.bg <- backgroundCorrect(array.raw, method='normexp', normexp.method='mle', offset=50)
		if (onecolor) {
			array <- array.bg$R
			array <- log2(array)
		} else {
			array.ma <- normalizeWithinArrays(array.bg, method=intra.array.normalization)
			array <- as.vector(array.ma$M)
			# take inverse for dye swaps
			if (cangem.samples[i, 'SampleChannel'] == 'Cy3')
				array <- -array
		}
		names(array) <- array.bg$genes[,cangem.samples[i,'ProbeNames']]
		# average replicate probes
		replicate.probes <- unique(names(array)[duplicated(names(array))])
		array.uniques <- array[!names(array) %in% replicate.probes]
		array.replicates <- array[names(array) %in% replicate.probes]
		array.replicates.avg <- aggregate(array.replicates, list(probe=names(array.replicates)), mean, na.rm=TRUE)
		array.replicates <- array.replicates.avg$x
		names(array.replicates) <- array.replicates.avg$probe
		array <- c(array.uniques, array.replicates)
		arrays[[i]] <- array
	}
	# go through all the loaded arrays and build a list of all probes found.
	all.probes <- character()
	for (i in 1:length(arrays))
		all.probes <- union(all.probes, names(arrays[[i]]))
	# build a matrix of the measurement values
	dat <- matrix(nrow=length(all.probes), ncol=length(chips), dimnames=list(all.probes, chips))
	for (i in 1:length(arrays))
		dat[,chips[i]] <- arrays[[i]][all.probes]
	# normalize between arrays
	dat <- normalizeBetweenArrays(dat, method=inter.array.normalization)
	dat <- round(dat, digits=2)
	dat <- as.data.frame(dat)
	chiptype <- cangem.samples$BioconductorPackage[1]
	if (is.na(chiptype)) {
		if (cangem.samples$Type[1] == 'miRNA') {
			chiptype <- 'miRNA'
		} else {
			chiptype <- 'cDNA'
		}
	}
} else if (cangem.samples$Format[1] == 'affymetrix') {
	library(affy)
	library(gcrma)
	onecolor <- TRUE
	for (i in 1:nrow(cangem.samples))
		download.file(paste('http://www.cangem.org/download.php?hybridization=', cangem.samples[i,'Accession'], auth, sep=''), cangem.samples[i,'FileName'], quiet=TRUE)
	raw <- ReadAffy(filenames=cangem.samples$FileName)
	chiptype<-paste(raw@annotation, '.db', sep='')
	calls<-as.data.frame(exprs(mas5calls(raw)))
	names(calls)<-paste("flag.", microarrays, sep="") # or names(calls) ???
	if (affymetrix.normalization == 'mas5')
		dat<-as.data.frame(round(log2(exprs(mas5(raw))), digits=2))
	rm(raw)
	gc()
	if (affymetrix.normalization == 'rma')
		dat <- as.data.frame(round(exprs(justRMA(filenames=cangem.samples$FileName)), digits=2))
	if (affymetrix.normalization == 'gcrma')
		dat <- as.data.frame(round(exprs(justGCRMA(filenames=cangem.samples$FileName, type='fullmodel', fast=TRUE, optimize.by='speed')), digits=2))
	names(dat) <- chips # or names(dat) ???
	dat <- data.frame(dat, calls)
	unlink(cangem.samples$FileName)
} else {
	stop('Unsupported file format: ', cangem.samples$Format[1])
}

phenodata <- data.frame(sample=microarrays)
phenodata$original_name <- cangem.samples$FileName
phenodata$chiptype <- chiptype
phenodata$group <- ''
phenodata$description <- cangem.samples$Name
cangem.samples$ProbeNames <- NULL
cangem.samples$BioconductorPackage <- NULL
if (onecolor) {
	cangem.samples$SampleChannel <- NULL
	cangem.samples$ReferenceSample <- NULL
	cangem.samples$ReferenceSex <- NULL
	cangem.samples$ReferenceAccession <- NULL
}

for (col in colnames(cangem.samples))
	if (all(is.na(cangem.samples[,col])))
		cangem.samples[,col] <- NULL

phenodata <- cbind(phenodata, cangem.samples)

dat2 <- data.frame(probe=rownames(dat), stringsAsFactors=FALSE)
rownames(dat2) <- dat2$probe

if (genome.build != 'none') {
	# load platform
	platform <- read.table(paste('http://www.cangem.org/download.php?platform=', cangem.samples$PlatformAccession[1], '&flag=', genome.build, auth, sep=''), sep='\t', header=TRUE, as.is=TRUE)
	colnames(platform) <- tolower(colnames(platform))
	colnames(platform)[colnames(platform)=='chr'] <- 'chromosome'
	rownames(platform) <- platform[,1]
	platform$chromosome <- factor(platform$chromosome, levels=c(1:22, "X", "Y", "MT"), ordered=TRUE)
	dat2 <- cbind(dat2, platform[dat2$probe, c('chromosome', 'start', 'end')])
	dat2 <- dat2[order(dat2$chromosome, dat2$start),]
}

if (chiptype != 'cDNA' && chiptype != 'miRNA') {
	# including gene names to data
	library(chiptype, character.only=T)
	lib2<-sub('.db','',chiptype)
	symbol<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep="")))))[dat2$probe,])
	genename<-gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep="")))))[dat2$probe,])
	symbol<-gsub("#", "", symbol)
	genename<-gsub("#", "", genename)
	dat2 <- cbind(dat2, symbol, description=genename)
}

dat2 <- cbind(dat2, dat[dat2$probe,])
if (ncol(dat)==1)
	colnames(dat2)[ncol(dat2)] <- chips[1]
dat2$probe <- NULL

write.table(dat2, file='normalized.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)
write.table(phenodata, file='phenodata.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)

# EOF
