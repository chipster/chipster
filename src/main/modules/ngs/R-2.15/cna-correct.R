# TOOL cna-correct.R: "Correct for GC content" (This tool corrects NGS data for GC content for copy number analysis.)
# INPUT cna-data-table.tsv: "Data table with log-transformed read counts" TYPE GENE_EXPRS
# OUTPUT corrected-cna-data-table.tsv: "Data table with log-transformed read counts corrected for GC content"

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-12-18

# add span parameter + family !!!

file <- 'cna-data-table.tsv'
input <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, as.is=TRUE, check.names=FALSE)
binbp <- input$end[1] - input$start[1] + 1
bin.size <- round(binbp / 1000)

bins <- input[,-grep("^chip\\.", colnames(input))]
counts <- as.matrix(input[,grep("^chip\\.", colnames(input))])
if (ncol(counts) == 1)
  colnames(counts) <- grep("^chip\\.", colnames(input), value=TRUE)
counts <- 2^counts - 1
dat <- list(bins=bins, counts=counts)

# update to new naming and new scale !!!
genome.build <- 'GRCh37'
gc <- read.table(file.path(chipster.tools.path, 'MPScall', genome.build, paste('gc.', bin.size, 'kbp.txt.gz', sep='')), header=TRUE, sep='\t', as.is=TRUE)
mappability <- read.table(file.path(chipster.tools.path, 'MPScall', genome.build, paste('mappability.', bin.size, 'kbp.txt.gz', sep='')), header=TRUE, sep='\t', as.is=TRUE)
gc$gc <- gc$gc * 100
gc <- gc[gc$chromosome != 'MT',]
mappability <- mappability[mappability$chromosome != 'MT',]
# or should gc and mappability be added elsewhere ???

dat[['bins']]$gc <- gc$gc
dat[['bins']]$mappability <- mappability$mappability

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

dat <- correct(dat)
dat2 <- cbind(dat[['bins']], log2(dat[['corrected']]+1), stringsAsFactors=FALSE)

# write outputs
options(scipen=10)
write.table(dat2, 'corrected-cna-data-table.tsv', quote=FALSE, sep='\t', na='')

# EOF
