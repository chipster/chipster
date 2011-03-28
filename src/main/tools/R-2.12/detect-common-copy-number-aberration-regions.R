# ANALYSIS "aCGH"/"Identify common regions from called aCGH data" (Reduces dimensionality of called aCGH data by identifying common breakpoints.)
# INPUT GENE_EXPRS aberrations.tsv
# OUTPUT regions.tsv, regions.pdf, regions-frequencies.pdf
# PARAMETER max.info.loss DECIMAL DEFAULT 0.01 (Maximal information loss allowed.)
# PARAMETER image.width INTEGER FROM 200 TO 6400 DEFAULT 2400 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 6400 DEFAULT 2400 (Height of the plotted network image)

# detect-common-copy-number-aberration-regions.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-03-28

library(CGHcall)
library(CGHregions)
library(WECCA)

dat <- read.table('aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

dat$chromosome[dat$chromosome=='X'] <- '23'
dat$chromosome[dat$chromosome=='Y'] <- '24'
dat$chromosome[dat$chromosome=='MT'] <- '25'
dat$chromosome <- as.integer(dat$chromosome)

calls <- as.matrix(dat[,grep("^flag\\.", names(dat))])
logratios <- as.matrix(dat[,grep("^chip\\.", names(dat))])

if (length(grep("^segmented\\.", names(dat)))>0) { # input contains probabilities (is a CGHcall object)
  segmented <- as.matrix(dat[,grep("^segmented\\.", names(dat))])
  probloss <- as.matrix(dat[,grep("^probloss\\.", names(dat))])
  probnorm <- as.matrix(dat[,grep("^probnorm\\.", names(dat))])
  probgain <- as.matrix(dat[,grep("^probgain\\.", names(dat))])
  probamp <- as.matrix(dat[,grep("^probamp\\.", names(dat))])

  if (ncol(probamp)==0) {
    cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=logratios, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
  } else {
    cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=logratios, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain, probamp=probamp), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
  }
  regions <- regioning(cgh, max.info.loss)
  dat2 <- data.frame(regions[['ann']])
  hardcalls <- regions[['hardcalls']]
  softcalls <- TRUE
} else { # input contains only calls and logratios, no probabilities
  cgh <- data.frame(Probe=row.names(dat), Chromosome=dat$chromosome, Start=dat$start, End=dat$end, calls)
  regions <- CGHregions(cgh, max.info.loss)
  dat2 <- data.frame(regions@featureData@data)
  hardcalls <- assayDataElement(regions, 'regions')
  softcalls <- FALSE
}

colnames(dat2)[1:5] <- c('chromosome', 'start', 'end', 'num.probes', 'ave.dist')
dat2$ave.dist <- NULL

# end column contains starting positions of the last probes of a region
# change them to the end positions of those probes
if (nrow(dat2[dat2$chromosome %in% dat$chromosome & dat2$end %in% dat$start,]) >= nrow(dat2[dat2$chromosome %in% dat$chromosome & dat2$end %in% dat$end,]))
  for (row in rownames(dat2))
    dat2[row, 'end'] <- dat[dat$chromosome == dat2[row, 'chromosome'] & dat$start == dat2[row, 'end'], 'end'][1]

dat2$loss.freq <- round(mean(as.data.frame(t(hardcalls==-1))), digits=3)
dat2$gain.freq <- round(mean(as.data.frame(t(hardcalls==1))), digits=3)
if (2 %in% hardcalls)
  dat2$amp.freq <- round(mean(as.data.frame(t(hardcalls==2))), digits=3)

dat2 <- cbind(dat2, hardcalls)
colnames(dat2) <- sub('calls\\.', 'flag\\.', colnames(dat2))

region.medians <- hardcalls
for (row in rownames(region.medians)) {
    region.medians[row,] <- apply(logratios[dat$chromosome == dat2[row, "chromosome"] &
                                            dat$start >= dat2[row, "start"] &
                                            dat$start <= dat2[row, "end"],
                                           ], 2, median)
}
colnames(region.medians) <- sub('flag\\.', 'chip\\.', colnames(region.medians))
dat2 <- cbind(dat2, region.medians)

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

write.table(dat2, file='regions.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# pdf(file='regions.pdf', width=image.width/72, height=image.height/72)
pdf(file='regions.pdf')
plot(regions)
dev.off()

# pdf(file='regions-frequencies.pdf', width=image.width/72, height=image.height/72)
pdf(file='regions-frequencies.pdf')
frequencyPlot(regions)
dev.off()

# EOF
