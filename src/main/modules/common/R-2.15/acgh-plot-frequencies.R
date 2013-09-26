# TOOL acgh-plot-frequencies.R: "Plot copy number aberration frequencies" (Plot copy number aberration frequencies for different groups.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT group-frequencies.pdf: group-frequencies.pdf 
# PARAMETER column: column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column defining sample groups.)
# PARAMETER chromosomes: chromosomes TYPE STRING DEFAULT 0 (The numbers of the chromosomes to be plotted, separated by commas. 0 means all chromosomes. Ranges are also supported (e.g. 1,3,7-10\).)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2013-04-04

source(file.path(chipster.common.path, 'CGHcallPlus.R'))

# read input files
file <- 'aberrations.tsv'
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, as.is=TRUE, check.names=FALSE)
phenodata <- read.table("phenodata.tsv", header=TRUE, sep="\t", as.is=TRUE, check.names=FALSE)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')



# build a cghCall object from the data
dat$chromosome[dat$chromosome=='X'] <- 23
dat$chromosome[dat$chromosome=='Y'] <- 24
dat$chromosome[dat$chromosome=='MT'] <- 25
dat$chromosome <- as.integer(dat$chromosome)

dat <- dat[dat$chromosome %in% 1:24,]

dat <- dat[order(dat$chromosome, dat$start),]

calls <- as.matrix(dat[,grep("flag", names(dat))])
copynumber <- as.matrix(dat[,grep("chip", names(dat))])
segmented <- as.matrix(dat[,grep("segmented", names(dat))])
probloss <- as.matrix(dat[,grep("probloss", names(dat))])
probnorm <- as.matrix(dat[,grep("probnorm", names(dat))])
probgain <- as.matrix(dat[,grep("probgain", names(dat))])
probamp <- as.matrix(dat[,grep("probamp", names(dat))])

if (ncol(segmented) == 0) { # no segments
  cgh <- new('cghRaw', copynumber=copynumber, featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
} else if (ncol(calls) == 0) { # segments, but no calls
  cgh <- new('cghSeg', assayData=assayDataNew(copynumber=copynumber, segmented=segmented), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
} else if (ncol(probnorm) == 0) { # calls, but no probabilities
  probloss <- calls == -1
  probnorm <- calls ==  0
  probgain <- calls ==  1
  cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
} else if (ncol(probamp) == 0) { # probabilities, but not for amplifications
  cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
} else { # probabilities even for amplifications
  cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain, probamp=probamp), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
}
sampleNames(cgh) <- phenodata$description

remove <- rowSums(is.na(copynumber)) + rowSums(is.na(segmented))
cgh <- cgh[remove == 0,]

# parse chromosomes to be plotted
chromosomes <- gsub('X', '23', chromosomes, ignore.case=TRUE)
chromosomes <- gsub('Y', '24', chromosomes, ignore.case=TRUE)
chromosomes <- gsub('MT', '25', chromosomes, ignore.case=TRUE)
chromosomes <- gsub('[^0-9,-]', ',', chromosomes)
items <- strsplit(chromosomes, ',')[[1]]
chrs.to.plot <- integer()
for (item in items) {
  item <- item[item!='']
  if (length(item)==0) next
  range <- strsplit(item, '-')[[1]]
  range <- range[range!='']
  if (length(range)==0) next
  chrs.to.plot <- c(chrs.to.plot, seq(range[1], range[length(range)]))
}
chrs.to.plot <- unique(chrs.to.plot)
chrs.to.plot <- chrs.to.plot[chrs.to.plot %in% dat$chromosome]
if (length(chrs.to.plot)==0)
  chrs.to.plot <- 0

if (0 %in% chrs.to.plot) {
  cond <- TRUE
} else {
  cond <- chromosomes(cgh) %in% chrs.to.plot
}

# plot
pdf(file='group-frequencies.pdf', paper='a4r', width=0, height=0)
if (column == 'EMPTY') {
  frequencyPlot(cgh[cond,])
} else {
  for (group in sort(unique(phenodata[,column]))) {
    if (is.na(group)) {
      frequencyPlot(cgh[cond,is.na(phenodata[,column])], main=paste(column, ' = ', group, sep=''))
    } else {
      frequencyPlot(cgh[cond,!is.na(phenodata[,column]) & phenodata[,column] == group], main=paste(column, ' = ', group, sep=''))
    }
  }
}
dev.off()

# EOF
