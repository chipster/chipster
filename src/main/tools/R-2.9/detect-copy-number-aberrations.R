# ANALYSIS "aCGH tools (beta testing)"/"Detect copy number aberrations" (Call copy number aberrations from aCGH log ratios.)
# INPUT GENE_EXPRS normalized.tsv
# OUTPUT aberrations.tsv, aberration-summary.png
# PARAMETER chromosomes INTEGER DEFAULT 23 (Number of chromosomes. Usually 23 for sex-matched reference samples and 22 otherwise.)
# PARAMETER normalization [median, mode, none] DEFAULT median (Normalization method.)
# PARAMETER cn.states [3, 4] DEFAULT 3 (Whether to call loss/normal/gain or loss/normal/gain/amplification.)

# detect-copy-number-aberrations.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-03-11

library(CGHcall)

dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
dat <- data.frame(probe=rownames(dat), dat[,c('chromosome', 'start', 'end')], dat[,grep("chip", names(dat))], stringsAsFactors=FALSE)

dat$chromosome[dat$chromosome=='X'] <- 23
dat$chromosome[dat$chromosome=='Y'] <- 24
dat$chromosome[dat$chromosome=='MT'] <- 25
dat$chromosome <- as.integer(dat$chromosome)

cgh <- make_cghRaw(dat)
cgh <- preprocess(cgh, nchrom=chromosomes)
cgh <- normalize(cgh, method=normalization)
cgh <- segmentData(cgh)
cgh <- postsegnormalize(cgh)
cgh <- CGHcall(cgh, nclass=as.integer(cn.states))

dat2 <- data.frame(probe=rownames(assayDataElement(cgh, "calls")), cgh@featureData@data, stringsAsFactors=FALSE)
colnames(dat2) <- c('probe', 'chromosome', 'start', 'end')

calls <- assayDataElement(cgh, 'calls')
colnames(calls) <- sub('chip.', 'flag.', colnames(calls))
dat2 <- cbind(dat2, calls)

copynumber <- assayDataElement(cgh, 'copynumber')
dat2 <- cbind(dat2, copynumber)

segmented <- assayDataElement(cgh, 'segmented')
colnames(segmented) <- sub('chip.', 'segmented.', colnames(segmented))
dat2 <- cbind(dat2, segmented)

probloss <- assayDataElement(cgh, 'probloss')
colnames(probloss) <- sub('chip.', 'probloss.', colnames(probloss))
dat2 <- cbind(dat2, probloss)

probnorm <- assayDataElement(cgh, 'probnorm')
colnames(probnorm) <- sub('chip.', 'probnorm.', colnames(probnorm))
dat2 <- cbind(dat2, probnorm)

probgain <- assayDataElement(cgh, 'probgain')
colnames(probgain) <- sub('chip.', 'probgain.', colnames(probgain))
dat2 <- cbind(dat2, probgain)

if (cn.states=='4') {
  probamp <- assayDataElement(cgh, 'probamp')
  colnames(probamp) <- sub('chip.', 'probamp.', colnames(probamp))
  dat2 <- cbind(dat2, probamp)
}

dat2$loss.freq <- mean(as.data.frame(t(assayDataElement(cgh, "calls")==-1)))
dat2$gain.freq <- mean(as.data.frame(t(assayDataElement(cgh, "calls")==1)))
if (cn.states=='4')
  dat2$amp.freq <- mean(as.data.frame(t(assayDataElement(cgh, "calls")==2)))

dat2$chromosome <- as.character(dat2$chromosome)
dat2$chromosome[dat2$chromosome=='23'] <- 'X'
dat2$chromosome[dat2$chromosome=='24'] <- 'Y'
dat2$chromosome[dat2$chromosome=='25'] <- 'MT'

write.table(dat2, file='aberrations.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)

bitmap(file='aberration-summary.png', width=600/72, height=600/72)
plot.summary(cgh)
dev.off()

# EOF