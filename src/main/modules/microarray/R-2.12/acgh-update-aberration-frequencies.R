# TOOL acgh-update-aberration-frequencies.R: "Update aberration frequencies for called copy number data" (Calculate frequencies of losses, gains and amplifications for called aCGH data. They are automatically calculated by the Call copy number aberrations from aCGH data tool, so normally there is no need to run this tool separately. But if you use e.g. the Extract samples from dataset tool, the frequencies should be updated using this script.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# OUTPUT aberration-frequencies.tsv: aberration-frequencies.tsv 

# calculate-aberration-frequencies.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-12-13

# read data set
dat <- read.table('aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
calls <- dat[,grep('flag', colnames(dat))]

# calculate aberration frequencies
dat$loss.freq <- round(rowMeans(calls == -1), digits=3)
dat$gain.freq <- round(rowMeans(calls == 1), digits=3)
if (2 %in% calls) {
  dat$amp.freq <- round(rowMeans(calls == 2), digits=3)
} else {
  dat$amp.freq <- NULL
}

write.table(dat, file='aberration-frequencies.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF
