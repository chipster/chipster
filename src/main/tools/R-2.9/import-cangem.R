# ANALYSIS Utilities/"Import from CanGEM" (Load a microarray data set from the CanGEM database)
# OUTPUT normalized.tsv, phenodata.tsv
# PARAMETER dataset STRING (data set accession)
# PARAMETER username STRING (username)
# PARAMETER password STRING (password)
# PARAMETER session STRING (session ID, can be used for authentication instead of a username/password combo)
# PARAMETER filtering [yes, no] DEFAULT yes (filter outliers)
# PARAMETER intra.array.normalization [none, median, loess] DEFAULT loess (normalization method)
# PARAMETER inter.array.normalization [none, quantile, scale] DEFAULT none (normalization method)
# PARAMETER map.to.build [GRCh37, NCBI36, NCBI35, none] DEFAULT GRCh37 (The genome build to use for adding the chromosome names and start and end base pair positions for the probes. Unless this is set to 'none', all probes that do not have a mapped position stored in CanGEM will be removed from the data.)

# import.from.cangem.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-02-17

library(limma)

# define filtering function
if (filtering == 'yes') {
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

# construct the string used in authenticating
if (session != '') {
  auth <- paste('&PHPSESSID=', session, sep='')
} else if (username != '' && password != '') {
  auth <- paste('&username=', username, '&password=', password, sep = '')
} else auth <- ''

# fetch list of arrays from CanGEM
samples <- read.table(paste('http://www.cangem.org/scripts/listhybs.php?dataset=', dataset, auth, sep=''),
  sep='\t', header=TRUE, as.is=TRUE)

# check for single platform, array type and file format
if (length(unique(samples$PlatformAccession)) > 1)
  stop('Multiple platforms used: ', paste(unique(samples$Platform), collapse=', '))
if (length(unique(samples$Type)) > 1)
  stop('Multiple array types used: ', paste(unique(samples$Type), collapse=', '))
if (length(unique(samples$Format)) > 1)
  stop('Multiple file formats used: ', paste(unique(samples$Format), collapse=', '))

arrays <- list()
names <- character()
onecolor <- FALSE
if (samples$Format[1] == 'agilent') {
  # load arrays one by one to be able to handle with files with different number of rows
  for (i in 1:nrow(samples)) {
    prob <- TRUE
    if (!onecolor) {
      # first try to read the 'normal' columns (means), but if it fails, try medians
      try({
        array.raw <- read.maimages(paste('http://www.cangem.org/download.php?hybridization=', samples[i,'Accession'], auth, sep=''),
          source=samples[i,'Format'], wt.fun=cangem.filter)
        prob <- FALSE
      }, silent=TRUE)
      if (prob) {
        try({
          array.raw <- read.maimages(paste('http://www.cangem.org/download.php?hybridization=', samples[i,'Accession'], auth, sep=''),
            source=samples[i,'Format'], columns=list(G='gMedianSignal', Gb='gBGMedianSignal', R='rMedianSignal', Rb='rBGMedianSignal'), wt.fun=cangem.filter)
          prob <- FALSE
        }, silent=TRUE)
      }
    }
    if (prob) {
      try({
        array.raw <- read.maimages(paste('http://www.cangem.org/download.php?hybridization=', samples[i,'Accession'], auth, sep=''),
          source=samples[i,'Format'], columns=list(G='gMedianSignal', Gb='gBGMedianSignal', R='gMedianSignal', Rb='gBGMedianSignal'), wt.fun=cangem.filter)
        prob <- FALSE
        onecolor <- TRUE
      }, silent=TRUE)
    }
    if (prob) {
       stop('Could not read file ', samples[i, 'Name'])
    }
    array.bg <- backgroundCorrect(array.raw, method='normexp', normexp.method='mle', offset=50)
    if (onecolor) {
      array <- array.bg$R
      array <- log2(array)
    } else {
      array.ma <- normalizeWithinArrays(array.bg, method=intra.array.normalization)
      array <- as.vector(array.ma$M)
      # take inverse for dye swaps
      if (samples[i, 'SampleChannel'] == 'Cy3')
        array <- -array
    }
    names(array) <- array.bg$genes[,samples[i,'ProbeNames']]
    # average replicate probes
    replicate.probes <- unique(names(array)[duplicated(names(array))])
    array.uniques <- array[!names(array) %in% replicate.probes]
    array.replicates <- array[names(array) %in% replicate.probes]
    array.replicates.avg <- aggregate(array.replicates, list(Probe=names(array.replicates)), mean, na.rm=TRUE)
    array.replicates <- array.replicates.avg$x
    names(array.replicates) <- array.replicates.avg$Probe
    array <- c(array.uniques, array.replicates)
    # dat[,samples[i, 'Name']] <- array[dat$Probe]
    arrays[[i]] <- array
    names[i] <- sprintf('microarray%.3i', i)
  }
} else {
  stop('Unsupported file format: ', samples$Format[1])
}

chips <- paste('chip.', names, sep='')

# go through all the loaded arrays and build a list of all probes found.
all.probes <- character()
for (i in 1:length(arrays))
  all.probes <- union(all.probes, names(arrays[[i]]))

# build a matrix of the measurement values
dat <- matrix(nrow=length(all.probes), ncol=length(chips), dimnames=list(all.probes, chips))
for (i in 1:length(arrays))
  dat[,chips[i]] <- arrays[[i]][all.probes]

dat <- normalizeBetweenArrays(dat, method=inter.array.normalization)
dat <- round(dat, digits=2)

if (map.to.build=='none') {
  dat <- data.frame(Probe=rownames(dat), dat)
} else {
  # load platform
  platform <- read.table(paste('http://www.cangem.org/download.php?platform=', samples$PlatformAccession[1], '&flag', map.to.build, auth, sep=''),
    sep='\t', header=TRUE, as.is=TRUE)
  rownames(platform) <- platform$Probe
  platform$Chr <- factor(platform$Chr, levels=c(1:22, "X", "Y", "MT"), ordered=TRUE)
  platform <- platform[order(platform$Chr, platform$Start),]
  dat <- data.frame(platform, dat[platform$Probe,])
}

if (onecolor) {
  samples$SampleChannel <- NULL
  samples$ReferenceSample <- NULL
  samples$ReferenceSex <- NULL
  samples$ReferenceAccession <- NULL
}
phenodata <- data.frame(sample=names)
phenodata$original_name <- samples$FileName
phenodata$chiptype = samples$BioconductorPackage
if (samples$Type[1] == 'miRNA Expression')
  phenodata$chiptype[is.na(phenodata$chiptype)] <- 'miRNA'
phenodata$chiptype[is.na(phenodata$chiptype)] <- 'cDNA'
phenodata$group = ''
phenodata$description <- samples$Name
samples$ProbeNames <- NULL
samples$BioconductorPackage <- NULL
phenodata <- cbind(phenodata, samples)

write.table(dat, file='normalized.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
write.table(phenodata, file='phenodata.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)

# EOF