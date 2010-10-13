# ANALYSIS "aCGH tools"/"Plot copy-number-induced gene expression" (Plot the expression levels of individual genes for a copy number vs. expression comparison. This tool must be run on the output from the tool Test for copy number induced expression changes - cn-induced-expression.tsv.)
# INPUT GENE_EXPRS cn-induced-expression.tsv
# OUTPUT cn-induced-expression-plot.png
# PARAMETER gene.ids STRING DEFAULT 1 (The gene.ids of the genes to be plotted, separated by commas. Ranges are also supported (e.g. 1,3,7-10).)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# plot-cn-induced-gene-expression.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2010-10-05

library(intCNGEan)

# read input file
dat <- read.table('cn-induced-expression.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

# if testing was done with analysis.type='regional', the resulting table contains four more columns
# those are removed so that we can calculate the number of samples from the number of columns
dat$reg.id <- NULL
dat$begin.reg <- NULL
dat$end.reg <- NULL
dat$shrinkage <- NULL

# return the data to the original order
# not sure if necessary, but just in case
dat <- dat[order(dat$gene.id),]
nosamp <- (ncol(dat)-12)/3
tuned <- list(datafortest=as.matrix(dat[,13:ncol(dat)]), lossorgain=dat$comparison, genestotest=dat$gene.id, callprobs=as.matrix(dat[,c('av.probs.1','av.probs.2')]), alleffects=dat$effect.size, ann=dat[,c('chromosome','start','end')], nosamp=nosamp)
rownames(tuned$datafortest) <- dat$probes
colnames(tuned$datafortest) <- sub('chip\\.', '', colnames(tuned$datafortest))
colnames(tuned$datafortest)[1:(2*nosamp)] <- ''
rownames(tuned$callprobs) <- dat$probes
colnames(tuned$callprobs) <- NULL
colnames(tuned$ann) <- c('Chromosome','Start','End')
names(tuned$lossorgain) <- dat$probes
names(tuned$genestotest) <- dat$probes

# parse the input string
gene.ids <- gsub('[^0-9,-]', ',', gene.ids)
items <- strsplit(gene.ids, ',')[[1]]
to.plot <- integer()
for (item in items) {
  item <- item[item!='']
  if (length(item)==0) next
  range <- strsplit(item, '-')[[1]]
  range <- range[range!='']
  if (length(range)==0) next
  to.plot <- c(to.plot, seq(range[1], range[length(range)]))
}
to.plot <- unique(to.plot)

# remove genes that are not present
to.plot <- intersect(to.plot, tuned$genestotest)

# check that we have something to plot
if (length(to.plot)==0)
  stop('CHIPSTER-NOTE: Nothing to plot.')

# plot
bitmap(file='cn-induced-expression-plot.png', width=image.width/72, height=image.height/72)
if (length(to.plot)==1) {
  intCNGEan.plot(gene.id=to.plot, tuned)
} else {
  sq <- sqrt(length(to.plot))
  rows <- ceiling(sq)
  cols <- ceiling(length(to.plot)/rows)
  par(mfrow=c(rows,cols))
  for (gene in to.plot)
    intCNGEan.plot(gene.id=gene, tuned)
}
dev.off()

# EOF