# ANALYSIS "aCGH"/"Call copy number aberrations from aCGH data" (Call copy number aberrations from aCGH log ratios.)
# INPUT GENE_EXPRS normalized.tsv
# OUTPUT aberrations.tsv, aberrations.png
# PARAMETER normalization [median, mode, none] DEFAULT none (Normalization method.)
# PARAMETER number.of.chromosomes INTEGER DEFAULT 23 (Number of chromosomes. Usually 23 for sex-matched reference samples and 22 otherwise.)
# PARAMETER number.of.copy.number.states [3, 4] DEFAULT 3 (Whether to call loss vs. normal vs. gain or loss vs. normal vs. gain vs. amplification.)
# PARAMETER minimum.number.of.probes.per.segment [2, 3, 4, 5] DEFAULT 2 (Minimum number of probes per segment.)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# detect-copy-number-aberrations.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2010-10-12

library(CGHcall)

dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

dat2 <- data.frame(probe=rownames(dat), dat[,c('chromosome', 'start', 'end')], dat[,grep("chip", names(dat))], stringsAsFactors=FALSE)
chips <- colnames(dat)[grep("chip", names(dat))]

if (ncol(dat2)==4)
  stop('CHIPSTER-NOTE: No array data found. The input file must have columns labeled with "chip.".')
if (ncol(dat2)==5)
  colnames(dat2)[5] <- chips[1]

dat2$chromosome[dat2$chromosome=='X'] <- 23
dat2$chromosome[dat2$chromosome=='Y'] <- 24
dat2$chromosome[dat2$chromosome=='MT'] <- 25
dat2$chromosome <- as.integer(dat2$chromosome)

cgh <- make_cghRaw(dat2)
cgh <- preprocess(cgh, nchrom=number.of.chromosomes)
cgh <- normalize(cgh, method=normalization)
cgh <- segmentData(cgh, min.width=as.integer(minimum.number.of.probes.per.segment))
cgh <- postsegnormalize(cgh)
cgh <- CGHcall(cgh, nclass=as.integer(number.of.copy.number.states))

dat3 <- data.frame(cgh@featureData@data)
colnames(dat3) <- c('chromosome', 'start', 'end')

for (col in c('cytoband', 'symbol', 'description', 'cnvs'))
  if (col %in% colnames(dat))
    dat3[,col] <- dat[rownames(dat3), col]

dat3$loss.freq <- round(mean(as.data.frame(t(assayDataElement(cgh, "calls")==-1))), digits=3)
dat3$gain.freq <- round(mean(as.data.frame(t(assayDataElement(cgh, "calls")==1))), digits=3)
if (number.of.copy.number.states=='4' && 2 %in% assayDataElement(cgh, 'calls'))
  dat3$amp.freq <- round(mean(as.data.frame(t(assayDataElement(cgh, "calls")==2))), digits=3)

calls <- assayDataElement(cgh, 'calls')
colnames(calls) <- sub('chip.', 'flag.', chips)
dat3 <- cbind(dat3, calls)

copynumber <- assayDataElement(cgh, 'copynumber')
colnames(copynumber) <- chips
dat3 <- cbind(dat3, copynumber)

segmented <- assayDataElement(cgh, 'segmented')
colnames(segmented) <- sub('chip.', 'segmented.', chips)
dat3 <- cbind(dat3, segmented)

probloss <- assayDataElement(cgh, 'probloss')
colnames(probloss) <- sub('chip.', 'probloss.', chips)
dat3 <- cbind(dat3, probloss)

probnorm <- assayDataElement(cgh, 'probnorm')
colnames(probnorm) <- sub('chip.', 'probnorm.', chips)
dat3 <- cbind(dat3, probnorm)

probgain <- assayDataElement(cgh, 'probgain')
colnames(probgain) <- sub('chip.', 'probgain.', chips)
dat3 <- cbind(dat3, probgain)

if (number.of.copy.number.states=='4') {
  probamp <- assayDataElement(cgh, 'probamp')
  colnames(probamp) <- sub('chip.', 'probamp.', chips)
  dat3 <- cbind(dat3, probamp)
}

dat3$chromosome <- as.character(dat3$chromosome)
dat3$chromosome[dat3$chromosome=='23'] <- 'X'
dat3$chromosome[dat3$chromosome=='24'] <- 'Y'
dat3$chromosome[dat3$chromosome=='25'] <- 'MT'

write.table(dat3, file='aberrations.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

bitmap(file='aberrations.png', width=image.width/72, height=image.height/72)
plot.summary(cgh)
dev.off()

# EOF