# TOOL acgh-plot-survival.R: "Plot survival curves for called copy number data" (Plotting of survival curves for called copy number data.)
# INPUT survival-test.tsv: survival-test.tsv TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT survival.pdf: survival.pdf
# PARAMETER survival: survival TYPE METACOLUMN_SEL DEFAULT survival (Phenodata column with survival data)
# PARAMETER status: status TYPE METACOLUMN_SEL DEFAULT status (Phenodata column with patient status: alive=0, dead=1)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-03-02

library(survival)

dat <- read.table('survival-test.tsv', header=TRUE, sep='\t', quote='', as.is=TRUE, row.names=1)
phenodata <- read.table('phenodata.tsv', header=TRUE, sep='\t')

s <- Surv(phenodata[,survival], phenodata[,status])
reg <- as.matrix(dat[, grep('^flag\\.', colnames(dat))])
call.legend <- c('loss', 'normal', 'gain')
names(call.legend) <- -1:1
call.cols <- c('red', 'black', 'blue')
names(call.cols) <- -1:1
pdf('survival.pdf', paper='a4r')
for (i in rownames(dat)) {
  f <- survfit(s ~ reg[i,])
  call <- sub('^.*=(.*)', '\\1', names(f$strata))
  n <- length(call)
  legend <- c(paste(call.legend[call], ' (', f$n, ')', sep=''), 'survival', '95% confidence interval')
  cols <- c(call.cols[call], 'black', 'black')
  ltys <- c(rep(0, n), 1, 2)
  pchs <- c(rep(22, n), NA, NA)
  if (!is.null(dat$cytoband)) {
    main <- paste('Survival for', dat[i, 'cytoband'])
  } else {
    main <- paste('Survival for region', i)
  }
  plot(f, main=main, xlab='t', ylab=expression(hat(S)(t)), col=cols)
  lines(f, conf.int='only', lty=2, col=cols)
  legend('bottomleft', legend=legend, lty=ltys, pch=pchs, pt.bg=cols, pt.cex=2, inset=0.01)
  mtext(paste(dat[i, 'chromosome'], ':', format(dat[i, 'start'], big.mark=','), '-', format(dat[i, 'end'], big.mark=','), sep=''), side=1, line=-1, adj=0.99)
  if (!is.null(dat$fdr))
    mtext(paste('FDR = ', sprintf('%.4f', dat[i, 'fdr'])), side=1, line=-2, adj=0.99)
  if (!is.null(dat$pvalue))
    mtext(paste('p = ', sprintf('%.4f', dat[i, 'pvalue'])), side=1, line=-3, adj=0.99)
}
dev.off()

# EOF
