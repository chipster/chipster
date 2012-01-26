# TOOL cna-plot-profile.R: "Plot copy number profile" (This tool plots a copy number profile from binned hits.)
# INPUT binned-hits.tsv: "Input file with binned read counts" TYPE GENERIC
# OUTPUT cna-profile.pdf: "Copy number profile"
# PARAMETER normalization: "Normalization" TYPE [median: median, mode: mode, none: none] DEFAULT median (Normalization method.)
# PARAMETER counts: "Counts" TYPE [count: "Original raw counts", corrected: "GC corrected counts"] DEFAULT corrected (Whether to use original raw read counts, or GC corrected ones.)
# PARAMETER chromosomes: chromosomes TYPE STRING DEFAULT 0 (The numbers of the chromosomes to be plotted, separated by commas. 0 means all chromosomes. Ranges are also supported (e.g. 1,3,7-10\).)
# PARAMETER resolution: resolution TYPE DECIMAL FROM 0.001 TO 1 DEFAULT 1 (Proportion of log-ratio data points to draw. Lower values lead to smaller file sizes and faster processing.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-12-22

source(file.path(chipster.tools.path, 'MPScall', 'CGHcallPlus-R-2.12.R'))

if (counts == 'count') {
  colClasses <- c('character', 'character', 'integer', 'integer', 'numeric', 'NULL')
} else {
  colClasses <- c('character', 'character', 'integer', 'integer', 'NULL', 'numeric')
}

dat <- read.table('binned-hits.tsv', header=TRUE, sep='\t', row.names=1, colClasses=colClasses)
dat[,4] <- dat[,4] - min(dat[,4], na.rm=TRUE) + 1 # to prevent negative values
dat[,4] <- log2(dat[,4])
dat <- data.frame(bin=1:nrow(dat), dat)

dat$chromosome[dat$chromosome=='X'] <- '23'
dat$chromosome[dat$chromosome=='Y'] <- '24'
dat$chromosome[dat$chromosome=='MT'] <- '25'
dat$chromosome <- as.integer(dat$chromosome)

cgh.raw <- make_cghRaw(dat)
cgh.pre <- preprocess(cgh.raw, nchrom=max(chromosomes(cgh.raw)))
cgh.nor <- normalize(cgh.pre, method=normalization)

sampleNames(cgh.nor) <- paste('CNA Profile With', nrow(dat), 'Bins')

# parse chromosomes to be plotted
chromosomes <- gsub('X', '23', chromosomes, ignore.case=TRUE)
chromosomes <- gsub('Y', '24', chromosomes, ignore.case=TRUE)
chromosomes <- gsub('MT', '25', chromosomes, ignore.case=TRUE)
chromosomes <- gsub('[^0-9,-]', ',', chromosomes)
items <- strsplit(chromosomes, ',')[[1]]
chrs.to.plot <- integer()
for (item in items) {
  item <- item[item!='']
  if (length(item)==0) next
  range <- strsplit(item, '-')[[1]]
  range <- range[range!='']
  if (length(range)==0) next
  chrs.to.plot <- c(chrs.to.plot, seq(range[1], range[length(range)]))
}
chrs.to.plot <- unique(chrs.to.plot)
chrs.to.plot <- chrs.to.plot[chrs.to.plot %in% dat$chromosome]
if (length(chrs.to.plot)==0)
  chrs.to.plot <- 0

# plot
pdf(file='cna-profile.pdf', paper='a4r', width=0, height=0)
if (0 %in% chrs.to.plot) {
  plot(cgh.nor[,1], dotres=1/resolution)
} else {
  plot(cgh.nor[chromosomes(cgh.nor) %in% chrs.to.plot, 1], dotres=1/resolution)
}
dev.off()

# EOF
