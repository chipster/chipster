# TOOL acgh-identify-regions.R: "Identify common regions from called copy number data" (Reduces dimensionality of called copy number data by identifying common breakpoints.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# OUTPUT regions.tsv: regions.tsv 
# OUTPUT region-frequencies.pdf: region-frequencies.pdf 
# PARAMETER max.info.loss: "maximum loss of information" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.01 (Maximal information loss allowed.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-01-11

source(file.path(chipster.common.path, 'CGHcallPlus.R'))
library(CGHregions)
library(WECCA)

file <- 'aberrations.tsv'
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, as.is=TRUE, check.names=FALSE)

dat$chromosome[dat$chromosome=='X'] <- '23'
dat$chromosome[dat$chromosome=='Y'] <- '24'
dat$chromosome[dat$chromosome=='MT'] <- '25'
dat$chromosome <- as.integer(dat$chromosome)

calls <- as.matrix(dat[,grep("^flag\\.", names(dat))])
logratios <- as.matrix(dat[,grep("^chip\\.", names(dat))])
segmented <- as.matrix(dat[,grep("^segmented\\.", names(dat))])

if (length(grep("^probnorm\\.", names(dat)))>0) { # input contains probabilities (is a CGHcall object)
  probdel <- as.matrix(dat[,grep("^probdel\\.", names(dat))])
  probloss <- as.matrix(dat[,grep("^probloss\\.", names(dat))])
  probnorm <- as.matrix(dat[,grep("^probnorm\\.", names(dat))])
  probgain <- as.matrix(dat[,grep("^probgain\\.", names(dat))])
  probamp <- as.matrix(dat[,grep("^probamp\\.", names(dat))])

  if (ncol(probamp)==0) { # CGHcall run with nclass=3
    cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=logratios, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
  } else if (ncol(probdel)==0) { # CGHcall run with nclass=4
    cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=logratios, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain, probamp=probamp), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
  } else { # CGHcall run with nclass=5
    cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=logratios, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain, probamp=probamp, probdloss=probdel), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
  }

  rm(calls, segmented, probdel, probloss, probnorm, probgain, probamp)
  gc()

  regions <- regioning(cgh, max.info.loss)
  dat2 <- data.frame(regions[['ann']])
  hardcalls <- regions[['hardcalls']]
  softcalls <- TRUE
} else { # input contains only calls and logratios, no probabilities
  cgh <- data.frame(Probe=row.names(dat), Chromosome=dat$chromosome, Start=dat$start, End=dat$end, calls)

  rm(calls, segmented)
  gc()

  regions <- CGHregions(cgh, max.info.loss)
  dat2 <- fData(regions)
  hardcalls <- regions(regions)
  softcalls <- FALSE
}

colnames(dat2)[1:5] <- c('chromosome', 'start', 'end', 'num.probes', 'ave.dist')
dat2$ave.dist <- NULL

# end column contains starting positions of the last probes/bin of a region
# change them to the end positions of those probes/bins
if (nrow(dat2[dat2$chromosome %in% dat$chromosome & dat2$end %in% dat$start,]) >= nrow(dat2[dat2$chromosome %in% dat$chromosome & dat2$end %in% dat$end,]))
  for (row in rownames(dat2))
    dat2[row, 'end'] <- dat[dat$chromosome == dat2[row, 'chromosome'] & dat$start == dat2[row, 'end'], 'end'][1]

# calculate frequencies
if (-2 %in% hardcalls)
  dat2$del.freq <- round(rowMeans(hardcalls == -2), digits=3)
dat2$loss.freq <- round(rowMeans(hardcalls == -1), digits=3)
dat2$gain.freq <- round(rowMeans(hardcalls == 1), digits=3)
if (2 %in% hardcalls)
  dat2$amp.freq <- round(rowMeans(hardcalls == 2), digits=3)

dat2 <- cbind(dat2, hardcalls)
colnames(dat2) <- sub('^calls\\.', 'flag.', colnames(dat2))

# calculate median logratios and segments
region.medians <- segmented.medians <- hardcalls
for (row in rownames(region.medians)) {
  index <- dat$chromosome == dat2[row, "chromosome"] &
           dat$start >= dat2[row, "start"] &
           dat$start <= dat2[row, "end"]
  region.medians[row,] <- apply(logratios[index,], 2, median, na.rm=TRUE)
  segmented.medians[row,] <- apply(logratios[index,], 2, median, na.rm=TRUE)
}
colnames(region.medians) <- sub('^flag\\.', 'chip.', colnames(region.medians))
colnames(segmented.medians) <- sub('^flag\\.', 'segmented.', colnames(segmented.medians))
dat2 <- cbind(dat2, region.medians, segmented.medians)

# append soft calls if we have them
if (softcalls) {
  colnames(regions[['softcalls']]) <- sub('_flag', '', colnames(regions[['softcalls']]))
  dat2 <- cbind(dat2, regions[['softcalls']])
  # create the cghRegions object for plotting
  regions <- new('cghRegions', assayData=assayDataNew(regions=regions[['hardcalls']]), featureData=new('AnnotatedDataFrame', data=regions[['ann']]))
}

dat2$chromosome <- as.character(dat2$chromosome)
dat2$chromosome[dat2$chromosome=='23'] <- 'X'
dat2$chromosome[dat2$chromosome=='24'] <- 'Y'
dat2$chromosome[dat2$chromosome=='25'] <- 'MT'

rownames(dat2) <- sprintf('%s:%i-%i', dat2$chromosome, dat2$start, dat2$end)

options(scipen=10)
write.table(dat2, file='regions.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

pdf(file='region-frequencies.pdf', paper='a4r', width=0, height=0)
frequencyPlot(regions)
plot(regions)
dev.off()

# EOF
