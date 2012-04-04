# TOOL cna-plot-profile.R: "Plot copy number profile" (Plot copy number profiles of individual samples.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT cna-profile.pdf: "Copy number profile"
# PARAMETER samples: samples TYPE STRING DEFAULT 1 (The numbers of the samples to be plotted, separated by commas. Ranges are also supported (e.g. 1,3,7-10\).)
# PARAMETER chromosomes: chromosomes TYPE STRING DEFAULT 0 (The numbers of the chromosomes to be plotted, separated by commas. 0 means all chromosomes. Ranges are also supported (e.g. 1,3,7-10\).)
# PARAMETER resolution: resolution TYPE DECIMAL FROM 0 TO 1 DEFAULT 1 (Proportion of data points to draw. Lower values lead to smaller file sizes and faster processing.)
# PARAMETER ploidy: "phenodata column with ploidy of samples" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column with ploidy)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-03-20

source(file.path(chipster.tools.path, 'MPScall', 'CGHcallPlus-R-2.12.R'))

# read input files
dat <- read.table('aberrations.tsv', header=TRUE, sep='\t', row.names=1, as.is=TRUE)
phenodata <- read.table("phenodata.tsv", header=TRUE, sep="\t", as.is=TRUE)

if (ploidy == 'EMPTY') {
  ploidy <- rep(2, nrow(phenodata))
} else {
  ploidy <- as.numeric(phenodata[,ploidy])
  ploidy[is.na(ploidy)] <- 2
}

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

if (0 %in% samples.to.plot)
  samples.to.plot <- 1:nrow(phenodata)

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

dat$chromosome[dat$chromosome=='X'] <- 23
dat$chromosome[dat$chromosome=='Y'] <- 24
dat$chromosome[dat$chromosome=='MT'] <- 25
dat$chromosome <- as.integer(dat$chromosome)

if (! 0 %in% chrs.to.plot)
  dat <- dat[dat$chromosome %in% chrs.to.plot,]

calls <- as.matrix(dat[,grep("^flag\\.", names(dat))])
cn <- as.matrix(dat[,grep("^copynumber\\.", names(dat))])
ratios <- as.matrix(dat[,grep("^chip\\.", names(dat))])
segmented <- as.matrix(dat[,grep("^segmented\\.", names(dat))])
ratios <- 2^ratios
segmented <- 2^segmented

.getChromosomeLengths <- function(build) {
    build <- as.integer(gsub('[^0-9]', '', build))
    if (build == 34 || build == 16) {
       chromosome.lengths <- c(246127941, 243615958, 199344050, 191731959, 181034922, 170914576, 158545518, 146308819, 136372045, 135037215, 134482954, 132078379, 113042980, 105311216, 100256656, 90041932, 81860266, 76115139, 63811651, 63741868, 46976097, 49396972, 153692391, 50286555)
    } else if (build == 35 || build == 17) {
       chromosome.lengths <- c(245522847, 243018229, 199505740, 191411218, 180857866, 170975699, 158628139, 146274826, 138429268, 135413628, 134452384, 132449811, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49554710, 154824264, 57701691)
    } else if (build == 36 || build == 18) {
       chromosome.lengths <- c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
    } else {
       chromosome.lengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
    }
    names(chromosome.lengths) <- 1:24
    chromosome.lengths
}

plot.freec <- 
function (x, i, dotres=1, ploidy=2, ylab='copy number', build="GRCh37",... )
{
    ylimit <- c(0, 3*ploidy)
    chrom <- x$chromosome
    pos <- x$start
    uni.chrom <- unique(chrom)
    chrom.lengths <- .getChromosomeLengths(build)[as.character(uni.chrom)]
    chrom.ends <- integer()
    cumul <- 0
    for (j in uni.chrom) {
        pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
        cumul <- cumul + chrom.lengths[as.character(j)]
        chrom.ends <- c(chrom.ends, cumul)
    }
    names(chrom.ends) <- names(chrom.lengths)
    nclone <- nrow(x)
    whichtoplot <- seq(1,nclone,by=dotres)
    cols <- c('red', 'black', 'blue')
    names(cols) <- -1:1
    # for (i in colnames(cn)) {
        plot(pos[whichtoplot], ploidy*ratios[whichtoplot,i], cex=.1, main=phenodata$description[i], ylab=ylab, xlab="chromosomes", ylim=ylimit, xaxt="n", xaxs="i", col='dimgray')
        if (dotres != 1)
            mtext(paste('Plot resolution: ', 100/dotres, '%', sep=''), side=3, line=0)
        abline(h=ploidy)
        if (length(chrom.ends) > 1)
            for (j in names(chrom.ends)[-length(chrom.ends)])
                abline(v=chrom.ends[j], lty='dashed')
        ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
        axis(side=1,at=ax,labels=uni.chrom,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
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
  plot.freec(dat, sample, dotres=1/resolution, ploidy=ploidy[sample])
dev.off()

# EOF
