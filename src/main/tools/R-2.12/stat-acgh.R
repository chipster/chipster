# ANALYSIS "aCGH"/"Group tests for called aCGH data" (Statistical tests between two or more groups for called aCGH data. The testing is recommended to be performed after running the Identify common regions from called aCGH data tool.)
# INPUT GENE_EXPRS regions.tsv, GENERIC phenodata.tsv
# OUTPUT groups-test.tsv
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER test.statistic [Chi-square, Wilcoxon, KW] DEFAULT Chi-square (The test to use: either Chi-square, Wilcoxon, or Kruskal-Wallis.)
# PARAMETER number.of.permutations INTEGER DEFAULT 10000 (The number of permutations. At least 10000 recommended for final calculations.)

# stat-acgh.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2010-10-06

library(CGHtest)

dat <- read.table('regions.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
phenodata <- read.table('phenodata.tsv', header=TRUE, sep='\t')

groupnames <- unique(phenodata[,column])
groupnames <- groupnames[!is.na(groupnames)]
groupnames <- groupnames[groupnames!='']

first.data.col <- min(grep('chip', names(dat)), grep('flag', names(dat)))
data.info <- dat[,1:(first.data.col-1)]

calls <- as.matrix(dat[,grep('flag', colnames(dat))])

datacgh <- data.frame()
group.sizes <- integer()
for (group in groupnames) {
  group.samples <- which(phenodata[,column] == group & !is.na(phenodata[,column]))
  group.calls <- calls[,group.samples]
  if (nrow(datacgh)==0) {
    datacgh <- group.calls
  } else {
    datacgh <- cbind(datacgh, group.calls)
  }
  group.sizes <- c(group.sizes, length(group.samples))
  data.info[,paste('loss.freq.', group, sep='')] <- round(mean(as.data.frame(t(group.calls==-1))), digits=3)
  data.info[,paste('gain.freq.', group, sep='')] <- round(mean(as.data.frame(t(group.calls==1))), digits=3)
}

pvs <-  pvalstest(datacgh, data.info, teststat=test.statistic, group=group.sizes, groupnames=groupnames, niter=number.of.permutations)
fdrs <- fdrperm(pvs)

write.table(fdrs, file="groups-test.tsv", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

# EOF