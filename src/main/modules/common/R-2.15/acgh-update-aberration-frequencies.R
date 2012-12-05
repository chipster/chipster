# TOOL acgh-update-aberration-frequencies.R: "Update aberration frequencies for called copy number data" (Calculate frequencies of losses, gains and amplifications for called aCGH data. They are automatically calculated by the Call copy number aberrations from aCGH data tool, so normally there is no need to run this tool separately. But if you use e.g. the Extract samples from dataset tool, the frequencies should be updated using this script.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# OUTPUT aberration-frequencies.tsv: aberration-frequencies.tsv 

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-10-12

# read data set
file <- 'aberrations.tsv'
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, as.is=TRUE, check.names=FALSE)
calls <- dat[,grep('flag', colnames(dat))]

# calculate aberration frequencies
dat$loss.freq <- round(rowMeans(calls == -1), digits=3)
dat$gain.freq <- round(rowMeans(calls == 1), digits=3)
if (2 %in% calls) {
  dat$amp.freq <- round(rowMeans(calls == 2), digits=3)
} else {
  dat$amp.freq <- NULL
}

options(scipen=10)
write.table(dat, file='aberration-frequencies.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF
