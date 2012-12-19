# TOOL cna-bin-hits.R: "Count reads in fixed bins" (Divide the genome into non-overlapping bins of fixed size, calculate reads in each bin, and correct for GC content.)
# INPUT alignment.bam: "Input BAM file" TYPE GENERIC
# OUTPUT binned-hits.tsv: "Binned hits"
# PARAMETER organism: "organism" TYPE [human: human] DEFAULT human (Organism.)
# PARAMETER genome.build: "genome build" TYPE [GRCh37: GRCh37, NCBI36: NCBI36] DEFAULT GRCh37 (Genome build.)
# PARAMETER bin.size: "bin size" TYPE [1: "1 kbp", 5: "5 kbp", 10: "10 kbp", 15: "15 kbp", 30: "30 kbp", 50: "50 kbp", 100: "100 kbp"] DEFAULT 30 (Bin size.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-12-18

options(scipen=10)
bin.size <- as.integer(bin.size)
binbp <- bin.size * 1000

.get.genome <- function(organism='human', genome.build='GRCh37') {
  genome.build <- as.integer(gsub('[^0-9]', '', genome.build))
  chromosomes <- c(1:22, 'X', 'Y')
  genomes <- list(GRCh37=c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566), NCBI36=c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954))
  names(genomes[['GRCh37']]) <- chromosomes
  names(genomes[['NCBI36']]) <- chromosomes
  if (genome.build == 37 || genome.build == 19) {
    genome <- genomes[['GRCh37']]
  } else if (genome.build == 36 || genome.build == 18) {
    genome <- genomes[['NCBI36']]
  } else {
    stop('Unknown genome build.')
  }
  genome
}

.create.bins <- function(binbp, ...) {
  genome <- .get.genome(...)
  start <- end <- count <- integer()
  for (chromosome in names(genome)) {
    chromosome.size <- genome[chromosome]
    chromosome.starts <- seq(from=1, to=chromosome.size, by=binbp)
    chromosome.ends <- chromosome.starts + binbp - 1
    chromosome.ends[length(chromosome.ends)] <- chromosome.size
    chromosome.breaks <- c(chromosome.starts, chromosome.size)
    start <- c(start, chromosome.starts)
    end <- c(end, chromosome.ends)
  }
  bins <- data.frame(chromosome=rep(names(genome), ceiling(genome / binbp)), start, end, stringsAsFactors=FALSE)
}

genome <- .get.genome('human', genome.build)
bins <- .create.bins(binbp, 'human', genome.build)
bins$count <- 0

system(paste(file.path(chipster.tools.path, 'samtools', 'samtools'), ' view -F 0x0404 -q 1 alignment.bam | cut -f3,4 | tr -d chr > hits.txt', sep=''))

skip <- 0
chunk <- 100000000
while(1) {
  hits <- as.data.frame(scan('hits.txt', what=list(chromosome=character(), pos=integer()), sep='\t', nmax=chunk, skip=skip, quiet=TRUE), stringsAsFactors=FALSE)
  if (nrow(hits) == 0)
    break
  for (chromosome in unique(hits$chromosome)) {
    if (!chromosome %in% names(genome))
      next
    chromosome.size <- genome[chromosome]
    chromosome.starts <- seq(from=1, to=chromosome.size, by=binbp)
    chromosome.breaks <- c(chromosome.starts, chromosome.size)
    count <- hist(hits[hits$chromosome==chromosome, 'pos'], breaks=chromosome.breaks, plot=FALSE)$count
    bins[bins$chromosome==chromosome, 'count'] <- bins[bins$chromosome==chromosome, 'count'] + count
  }
  skip <- skip + chunk
  rm(hits, count); gc(FALSE)
}
unlink('hits.txt')

rownames(bins) <- paste(bins$chromosome, ':', bins$start, '-', bins$end, sep='')
write.table(bins, 'binned-hits.tsv', quote=FALSE, sep='\t', na='')

# EOF
