# ANALYSIS "aCGH tools (beta testing)"/"Detect common copy number aberration regions" (Reduces dimensionality of called aCGH data by identifying common breakpoints.)
# INPUT GENE_EXPRS aberrations.tsv
# OUTPUT regions.tsv, aberration-regions.png, aberration-frequencies.png
# PARAMETER max.info.loss DECIMAL DEFAULT 0.01 (Maximal information loss allowed.)

# detect-common-copy-number-aberration-regions.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-03-11

library(CGHcall)
library(CGHregions)

dat <- read.table('aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

dat$chromosome[dat$chromosome=='X'] <- 23
dat$chromosome[dat$chromosome=='Y'] <- 24
dat$chromosome[dat$chromosome=='MT'] <- 25
dat$chromosome <- as.integer(dat$chromosome)

calls <- as.matrix(dat[,grep("flag", names(dat))])
copynumber <- as.matrix(dat[,grep("chip", names(dat))])
segmented <- as.matrix(dat[,grep("segmented", names(dat))])
probloss <- as.matrix(dat[,grep("probloss", names(dat))])
probnorm <- as.matrix(dat[,grep("probnorm", names(dat))])
probgain <- as.matrix(dat[,grep("probgain", names(dat))])
probamp <- as.matrix(dat[,grep("probamp", names(dat))])

if (ncol(probamp)==0) {
  cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
} else {
  cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain, probamp=probamp), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))
}

regions <- CGHregions(cgh, max.info.loss)

dat2 <- data.frame(regions@featureData@data, assayDataElement(regions, 'regions'))
colnames(dat2)[1:5] <- c('chromosome', 'start', 'end', 'num.probes', 'ave.dist')
names(dat2) <- sub('calls.', 'flag.', names(dat2))

# end column contains starting positions of the last probes of a region
# change them to the end positions of those probes
if (nrow(dat2[dat2$chromosome %in% dat$chromosome & dat2$end %in% dat$start,]) >= nrow(dat2[dat2$chromosome %in% dat$chromosome & dat2$end %in% dat$end,]))
  for (row in rownames(dat2))
    dat2[row, 'end'] <- dat[dat$chromosome == dat2[row, 'chromosome'] & dat$start == dat2[row, 'end'], 'end'][1]

dat2$chromosome <- as.character(dat2$chromosome)
dat2$chromosome[dat2$chromosome=='23'] <- 'X'
dat2$chromosome[dat2$chromosome=='24'] <- 'Y'
dat2$chromosome[dat2$chromosome=='25'] <- 'MT'

dat2$loss.freq <- mean(as.data.frame(t(assayDataElement(regions, 'regions')==-1)))
dat2$gain.freq <- mean(as.data.frame(t(assayDataElement(regions, 'regions')==1)))
dat2$amp.freq <- mean(as.data.frame(t(assayDataElement(regions, 'regions')==2)))

write.table(dat2, file='regions.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

bitmap(file='aberration-regions.png', width=600/72, height=600/72)
plot(regions)
dev.off()

bitmap(file='aberration-frequencies.png', width=600/72, height=600/72)
frequencyPlot(regions)
dev.off()

# EOF
