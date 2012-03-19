# TOOL cna-bin-hits.R: "Count reads in fixed bins" (Divide the genome into non-overlapping bins of fixed size, calculate reads in each bin, and correct for GC content.)
# INPUT alignment.bam: "Input BAM file" TYPE GENERIC
# OUTPUT binned-hits.tsv: "Binned hits"
# PARAMETER organism: "Organism" TYPE [human: human] DEFAULT human (Organism.)
# PARAMETER genome.build: "Genome build" TYPE [GRCh37: GRCh37] DEFAULT GRCh37 (Genome build.)
# PARAMETER bin.size: "Bin size" TYPE [1: "1 kbp", 5: "5 kbp", 10: "10 kbp", 15: "15 kbp", 30: "30 kbp", 50: "50 kbp", 100: "100 kbp"] DEFAULT 30 (Bin size.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-03-01

library(limma)
bin.size <- as.integer(bin.size)
binbp <- bin.size * 1000
binstring <- paste(bin.size, 'kbp', sep='')

# chipster.tools.path <- "/home/ischeini/.opt/"
system(paste(file.path(chipster.tools.path, 'samtools', 'samtools'), ' view -F 0x0404 -q 1 alignment.bam | cut -f3,4 | tr -d chr > hits.txt', sep=''))

chromosomes <- c(1:22, 'X', 'Y', 'MT')
genomes <- list(GRCh37=c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 16571)) # Y not 59034049
names(genomes[['GRCh37']]) <- chromosomes
start <- end <- count <- integer()
for (chromosome in chromosomes) {
  chromosome.size <- genomes[['GRCh37']][chromosome]
  chromosome.starts <- seq(from=1, to=chromosome.size, by=binbp)
  chromosome.ends <- chromosome.starts + binbp - 1
  chromosome.ends[length(chromosome.ends)] <- chromosome.size
  chromosome.breaks <- c(chromosome.starts, chromosome.size)
  start <- c(start, chromosome.starts)
  end <- c(end, chromosome.ends)
}
bins <- data.frame(chromosome=rep(chromosomes, ceiling(genomes[['GRCh37']] / binbp)), start, end, count=0, stringsAsFactors=FALSE)
f <- file('hits.txt')
skip <- 0
chunk <- 100000000
while(1) {
  hits <- as.data.frame(scan('hits.txt', what=list(chromosome=character(), pos=integer()), sep='\t', nmax=chunk, skip=skip, quiet=TRUE))
  if (nrow(hits) == 0)
    break
  for (chromosome in unique(hits$chromosome)) {
    if (!chromosome %in% chromosomes)
      next
    chromosome.size <- genomes[['GRCh37']][chromosome]
    chromosome.starts <- seq(from=1, to=chromosome.size, by=binbp)
    chromosome.breaks <- c(chromosome.starts, chromosome.size)
    count <- hist(hits[hits$chromosome==chromosome, 'pos'], breaks=chromosome.breaks, plot=FALSE)$count
    bins[bins$chromosome==chromosome, 'count'] <- bins[bins$chromosome==chromosome, 'count'] + count
  }
  skip <- skip + chunk
  rm(hits, count); gc(FALSE)
}
close(f)
unlink('hits.txt') 

gc.correct <- function(x) {
  # x$count[x$count==0] <- NA # ?
  means <- numeric(length=length(gc.intervals))
  for (i in 1:length(gc.intervals))
    means[i] <- mean(x[!is.na(gc$permil) & gc$permil==gc.intervals[i],'count'], na.rm=TRUE)
  loess.means <- loessFit(means, gc.intervals)$fitted
  correction <- median(loess.means, na.rm=TRUE) - loess.means
  for (i in 1:length(gc.intervals))
    x[!is.na(gc$permil) & gc$permil==gc.intervals[i],'corrected'] <- x[!is.na(gc$permil) & gc$permil==gc.intervals[i],'count'] + correction[i]
  x$corrected <- round(x$corrected)
  # x$corrected - min(x$corrected, na.rm=TRUE) # to prevent negative values
}

# gc <- read.table(file.path(chipster.tools.path, '..', '..', '.MPScall', paste('gc.', bin.size, 'kbp.txt.gz', sep='')), header=TRUE, sep='\t', as.is=TRUE, colClasses=c('character', 'integer', 'integer', 'numeric'))
gc <- read.table(file.path(chipster.tools.path, 'MPScall', genome.build, paste('gc.', bin.size, 'kbp.txt.gz', sep='')), header=TRUE, sep='\t', as.is=TRUE, colClasses=c('character', 'integer', 'integer', 'numeric'))
gc$permil <- as.integer(round(gc$gc*1000))
gc.intervals <<- seq(from=range(gc$permil, finite=TRUE)[1], to=range(gc$permil, finite=TRUE)[2])
bins$corrected <- gc.correct(bins)

rownames(bins) <- paste('bin-', 1:nrow(bins), sep='')
options(scipen=10)
write.table(bins, 'binned-hits.tsv', quote=FALSE, sep='\t', na='')

# EOF
