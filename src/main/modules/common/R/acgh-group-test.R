# TOOL acgh-group-test.R: "Group tests for called copy number data" (Statistical tests between two or more groups for called copy number data. The testing is recommended to be performed after running the Identify common regions from called copy number data tool.)
# INPUT regions.tsv: regions.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT groups-test.tsv: groups-test.tsv 
# OUTPUT groups-test.pdf: groups-test.pdf 
# PARAMETER column: Column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER test.statistic: "Test statistic" TYPE [Chi-square: Chi-square, Wilcoxon: Wilcoxon, KW: KW] DEFAULT Chi-square (The test to use: either Chi-square, Wilcoxon, or Kruskal-Wallis.)
# PARAMETER number.of.permutations: "Number of permutations" TYPE INTEGER DEFAULT 10000 (The number of permutations. At least 10000 recommended for final calculations.)
# PARAMETER test.aberrations: "Test aberrations" TYPE [1: gains, -1: losses, 0: both] DEFAULT 0 (Whether to test only for gains or losses, or both.) 

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-23

source(file.path(chipster.common.path, 'library-Chipster.R'))
library(gtools)

dat <- readData("regions.tsv")
phenodata <- readPhenodata("phenodata.tsv")

groupnames <- unique(phenodata[,column])
groupnames <- groupnames[!is.na(groupnames)]
groupnames <- groupnames[groupnames!='']

data.info <- dat[, annotationColumns(dat)]
calls <- as.matrix(dat[,grep('^flag\\.', colnames(dat))])

datacgh <- data.frame()
freqs <- data.frame(row.names=rownames(dat))
group.sizes <- integer()
for (group in groupnames) {
  group.samples <- which(phenodata[,column] == group & !is.na(phenodata[,column]))
  group.calls <- calls[,group.samples, drop=FALSE]
  if (nrow(datacgh)==0) {
    datacgh <- group.calls
  } else {
    datacgh <- cbind(datacgh, group.calls)
  }
  group.sizes <- c(group.sizes, length(group.samples))
  if (-2 %in% calls)
    freqs[,paste('del.freq.', group, sep='')] <- round(rowMeans(group.calls == -2), digits=3)
  freqs[,paste('loss.freq.', group, sep='')] <- round(rowMeans(group.calls == -1), digits=3)
  freqs[,paste('gain.freq.', group, sep='')] <- round(rowMeans(group.calls == 1), digits=3)
  if (2 %in% calls)
    freqs[,paste('amp.freq.', group, sep='')] <- round(rowMeans(group.calls == 2), digits=3)
}

# first try parallel computing
prob <- TRUE
try({
  library(CGHtestpar)
  pvs <-  pvalstest(datacgh, data.info, teststat=test.statistic, group=group.sizes, groupnames=groupnames, lgonly=as.integer(test.aberrations), niter=number.of.permutations, ncpus=4)
  fdrs <- fdrperm(pvs)
  prob <- FALSE
}, silent=TRUE)
# if problems, fall back to sequential computing
if (prob) {
  library(CGHtest)
  pvs <-  pvalstest(datacgh, data.info, teststat=test.statistic, group=group.sizes, groupnames=groupnames, lgonly=as.integer(test.aberrations), niter=number.of.permutations)
  fdrs <- fdrperm(pvs)
}

fdrs <- cbind(fdrs, freqs, dat[, dataColumns(dat)])

writeData(fdrs, "groups-test.tsv")

FDRplot <- function(fdrs, which, a, b, main = 'Frequency Plot with FDR',...) {
  par(mar=c(5,4,4,5) + 0.1)
  cols <- c('blue', 'red')
  names(cols) <- c('gain', 'loss')
  chromosomes <- fdrs$chromosome
  a.freq <- fdrs[,paste0(which, '.freq.', a)]
  b.freq <- fdrs[,paste0(which, '.freq.', b)]
  fdr <- fdrs$fdr
  if ('num.probes' %in% colnames(fdrs)) {
    chromosomes <- rep(chromosomes, fdrs$num.probes)
    a.freq <- rep(a.freq, fdrs$num.probes)
    b.freq <- rep(b.freq, fdrs$num.probes)
    fdr <- rep(fdr, fdrs$num.probes)
  }
  plot(a.freq, ylim=c(-1,1), type='h', col=cols[which], xlab='chromosomes', ylab='frequency', xaxt='n', yaxt='n', main=main, ...)
  points(-b.freq, type='h', col=cols[which])
  abline(h=0)
  abline(v=0, lty='dashed')
  chr.num <- chromosomes
  chr.num[chr.num == 'X'] <- 23
  chr.num[chr.num == 'Y'] <- 24
  chr.num[chr.num == 'MT'] <- 25
  chr.num <- as.integer(chr.num)
  chr.ends <- cumsum(table(chr.num))
  for(i in chr.ends)
    abline(v=i, lty='dashed')
  ax <- (chr.ends + c(0, chr.ends[-length(chr.ends)])) / 2
  axis(side=1, at=ax, labels=unique(chromosomes))
  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
  logfdr <- -log10(fdr)
  logfdr[logfdr == Inf] <- 10
  points(logfdr - 1, type='l')
  labels <- c(0.01, 0.05, 0.025, 0.1, 0.25, 0.5, 1)
  axis(side=4, at=-log10(labels) - 1, labels=labels, las=1)
  mtext('FDR', side=4, line=3)
  mtext(a, side=2, line=3, at=0.5)
  mtext(b, side=2, line=3, at=-0.5)
}

combs <- combinations(n=length(groupnames), r=2, v=groupnames)
pdf('groups-test.pdf', paper='a4r', width=0, height=0)
if (test.aberrations != '-1')
  for (i in seq_len(nrow(combs)))
    FDRplot(fdrs, 'gain', combs[i, 1], combs[i, 2])
if (test.aberrations != '1')
  for (i in seq_len(nrow(combs)))
    FDRplot(fdrs, 'loss', combs[i, 1], combs[i, 2])
dev.off()

# EOF
