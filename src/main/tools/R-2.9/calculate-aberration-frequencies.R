# ANALYSIS "aCGH tools (beta testing)"/"Calculate aberration frequencies for called aCGH data" (Calculate frequencies of losses, gains and amplifications for called aCGH data. They are automatically calculated by the Call copy number aberrations from aCGH data module, so normally there is no need to run this module separately. But if you use e.g. the Extract samples from dataset module, the frequencies should be updated using this script.)
# INPUT GENE_EXPRS aberrations.tsv
# OUTPUT aberration-frequencies.tsv

# calculate-aberration-frequencies.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-08-19

# read data set
dat <- read.table('aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
calls <- dat[,grep('flag', colnames(dat))]

# calculate aberration frequencies
dat$loss.freq <- round(mean(as.data.frame(t(calls==-1))), digits=3)
dat$gain.freq <- round(mean(as.data.frame(t(calls==1))), digits=3)
if (2 %in% calls) {
  dat$amp.freq <- round(mean(as.data.frame(t(calls==2))), digits=3)
} else {
  dat$amp.freq <- NULL
}

write.table(dat, file='aberration-frequencies.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF