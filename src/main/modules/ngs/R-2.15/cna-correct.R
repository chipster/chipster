# TOOL cna-correct.R: "Correct for GC content" (This tool corrects NGS data for GC content for copy number analysis.)
# INPUT cna-data-table.tsv: "Data table with log-transformed read counts" TYPE GENE_EXPRS
# OUTPUT corrected-cna-data-table.tsv: "Data table with log-transformed read counts corrected for GC content"
# PARAMETER span: "degree of smoothing" TYPE DECIMAL DEFAULT 0.65 (The alpha/span parameter for the degree of smoothing, controls the neighboorhood used for the fitting. For alpha < 1, the neighbourhood includes proportion alpha of the points, and these have tricubic weighting. For alpha > 1, all points are used, with the ‘maximum distance’ assumed to be alpha^1/2.)
# PARAMETER family: family TYPE [gaussian: gaussian, symmetric: symmetric] DEFAULT gaussian (Family used for the fitting. For gaussian, fittins is done by weighted least squares. For symmetric, a few iterations of an M-estimation procedure with Tukey's biweight are used.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2013-03-20

file <- 'cna-data-table.tsv'
input <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, as.is=TRUE, check.names=FALSE)

bins <- input[,-grep("^chip\\.", colnames(input))]
counts <- as.matrix(input[,grep("^chip\\.", colnames(input))])
if (ncol(counts) == 1)
  colnames(counts) <- grep("^chip\\.", colnames(input), value=TRUE)
counts <- 2^counts - 1
dat <- list(bins=bins, counts=counts)

correct <- function(dat, family='gaussian', span=0.35, ...) {
  dat[['corrected']] <- dat[['counts']]
  dat[['corrected']][] <- NA
  gc <- round(dat[['bins']]$gc)
  median.counts <- aggregate(dat[['counts']][dat[['bins']]$chromosome %in% 1:22,], by=list(gc=gc[dat[['bins']]$chromosome %in% 1:22]), median)
  median.counts <- median.counts[!is.na(median.counts$gc),]
  rownames(median.counts) <- median.counts$gc
  residuals <- matrix(nrow=nrow(median.counts), ncol=ncol(dat[['counts']]), dimnames=list(rownames(median.counts), colnames(dat[['counts']])))
  for (i in 1:ncol(dat[['counts']])) {
    vals <- median.counts[,i+1]
    l <- loess(vals ~ median.counts$gc, family=family, span=span) #, ...)
    correction <- median(l$fitted, na.rm=TRUE) - l$fitted
    names(correction) <- median.counts$gc
    corvals <- dat[['counts']][,i] + correction[as.character(gc)]
    corvals <- corvals - min(corvals, na.rm=TRUE) # to prevent negative values
    dat[['corrected']][,i] <- corvals
    residuals[,i] <- l$residuals
  }
  dat[['bins']]$mean.residual <- apply(residuals, 1, mean)[as.character(gc)]
  dat
}

dat <- correct(dat, span=span, family=family)
dat2 <- cbind(dat[['bins']], log2(dat[['corrected']]+1), stringsAsFactors=FALSE)

# write outputs
options(scipen=10)
write.table(dat2, 'corrected-cna-data-table.tsv', quote=FALSE, sep='\t', na='')

# EOF
