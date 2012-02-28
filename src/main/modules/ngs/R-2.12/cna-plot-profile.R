# TOOL cna-plot-profile.R: "Plot copy number profile" (Plot copy number profiles of individual samples.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT cna-profile.pdf: "Copy number profile"
# PARAMETER samples: samples TYPE STRING DEFAULT 1 (The numbers of the samples to be plotted, separated by commas. Ranges are also supported (e.g. 1,3,7-10\).)
# PARAMETER chromosomes: chromosomes TYPE STRING DEFAULT 0 (The numbers of the chromosomes to be plotted, separated by commas. 0 means all chromosomes. Ranges are also supported (e.g. 1,3,7-10\).)
# PARAMETER resolution: resolution TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.1 (Proportion of log-ratio data points to draw. Lower values lead to smaller file sizes and faster processing.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-02-28

source(file.path(chipster.tools.path, 'MPScall', 'CGHcallPlus-R-2.12.R'))

# read input files
dat <- read.table('aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
phenodata <- read.table("phenodata.tsv", header=TRUE, sep="\t", as.is=TRUE)

# parse samples to be plotted
if (samples=='0')
  samples <- paste('1-', nrow(phenodata), sep='')
samples <- gsub('[^0-9,-]', ',', samples)
items <- strsplit(samples, ',')[[1]]
samples.to.plot <- integer()
for (item in items) {
  item <- item[item!='']
  if (length(item)==0) next
  range <- strsplit(item, '-')[[1]]
  range <- range[range!='']
  if (length(range)==0) next
  samples.to.plot <- c(samples.to.plot, seq(range[1], range[length(range)]))
}
samples.to.plot <- unique(samples.to.plot)

# remove samples that are out of bounds
samples.to.plot <- samples.to.plot[samples.to.plot<=nrow(phenodata)]

# check that we have something to plot
if (length(samples.to.plot)==0)
  stop('CHIPSTER-NOTE: Nothing to plot.')

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

if (! 0 %in% chrs.to.plot)
  dat <- dat[dat$chromosome %in% chrs.to.plot,]

dat$chromosome[dat$chromosome=='X'] <- 23
dat$chromosome[dat$chromosome=='Y'] <- 24
dat$chromosome[dat$chromosome=='MT'] <- 25
dat$chromosome <- as.integer(dat$chromosome)

calls <- as.matrix(dat[,grep("^flag\\.", names(dat))])
cn <- as.matrix(dat[,grep("^copynumber\\.", names(dat))])
ratios <- as.matrix(dat[,grep("^chip\\.", names(dat))])
segmented <- as.matrix(dat[,grep("^segmented\\.", names(dat))])

plot.freec <- 
function (x, i, dotres=1, ploidy=2, ylab='copy number', build="GRCh37",... )
{
    ylimit <- c(0, 3*ploidy)
    chrom <- x$chromosome
    pos <- x$start
    chrom.ends <- .getCumulativeChromosomeEnds(build)[1:max(chrom)]
    nclone <- nrow(x)
    whichtoplot <- seq(1,nclone,by=dotres)
    cols <- c('red', 'black', 'blue')
    names(cols) <- -1:1
    for (j in 2:max(chrom))
        pos[chrom == j] <- pos[chrom == j] + chrom.ends[j-1]
    # for (i in colnames(cn)) {
        plot(pos[whichtoplot], ploidy*ratios[whichtoplot,i], cex=.1, main=i, ylab=ylab, xlab="chromosomes", ylim=ylimit, xaxt="n", xaxs="i", col='dimgray')
        if (dotres != 1)
            mtext(paste('Plot resolution: 1/',dotres, sep=''), side=3, line=0)
        # abline(h=ploidy)
        for (j in 2:max(chrom))
            abline(v=chrom.ends[j-1], lty='dashed')
        ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
        axis(side=1,at=ax,labels=unique(chrom),cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
        points(pos[whichtoplot], cn[whichtoplot,i], cex=.3, col=cols[as.character(calls[whichtoplot,i])])
        amps <- cn[,i]
        amps[amps >  ylimit[2]] <- ylimit[2]+0.15
        amps[amps <= ylimit[2]] <- NA
        points(pos, amps, pch=24, col='blue', bg='blue', cex=0.5)
        dels <- cn[,i]
        dels[dels <  ylimit[1]] <- ylimit[1]-0.15
        dels[dels >= ylimit[1]] <- NA
        points(pos, dels, pch=25, col='red', bg='red', cex=0.5)
        ### MAD
        mad.value <- round(mad(ratios[chrom < 23,i], na.rm=TRUE), digits=2)
        mtext(paste('MAD =', mad.value), side=3, line=0, adj=1)
        ### number of data points
        str <- paste(round(nclone / 1000), 'k x ', sep='')
        probe <- median(x$end - x$start + 1)
        if (probe < 1000) {
            str <- paste(str, probe, ' bp', sep='')
        } else {
            str <- paste(str, round(probe / 1000), ' kbp', sep='')
        }
        mtext(str, side=3, line=0, adj=0)
    # }
}

# plot
pdf(file='cna-profile.pdf', paper='a4r', width=0, height=0)
for (sample in samples.to.plot)
  if (0 %in% chrs.to.plot) {
    plot(cgh[,sample], dotres=1/resolution)
  } else {
    plot(cgh[chromosomes(cgh) %in% chrs.to.plot, sample], dotres=1/resolution)
  }
dev.off()

# EOF
