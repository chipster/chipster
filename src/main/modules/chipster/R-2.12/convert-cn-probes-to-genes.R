# TOOL convert-cn-probes-to-genes.R: "Convert called aCGH data from probes to genes" (Using the chromosomal locations of the probes of an aCGH data set, convert the data from probe-based to gene-based. The probes from which the copy number call for a given gene are determined as follows. If there are probes overlapping with the position of the gene, they are used. In case of no overlaps, the last preceding and first tailing probe are used.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# OUTPUT gene-aberrations.tsv: gene-aberrations.tsv 
# PARAMETER method.for.calls: method.for.calls TYPE [majority: majority, unambiguous: unambiguous] DEFAULT majority (The method majority means that if more than 50% of these probes give an aberrated signal, that call is used for the gene. The unambiguous method requires that all of the probes have the same call, otherwise the gene will be labeled as normal.)
# PARAMETER method.for.others: method.for.others TYPE [mean: mean, median: median] DEFAULT mean (Whether to use the mean or the median for calculating other data than copy number calls.)
# PARAMETER genome.build: genome.build TYPE [GRCh37: GRCh37, NCBI36: NCBI36, NCBI35: NCBI35, NCBI34: NCBI34] DEFAULT GRCh37 (The genome build to use for fetching the gene coordinates.)

# convert-cn-probes-to-genes.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-02-15

dat <- read.table('aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

# load genes
# genes <- read.table(paste('http://www.cangem.org/download.php?platform=CG-PLM-26&flag=', genome.build, sep=''), sep='\t', header=TRUE, row.names=1, as.is=TRUE)
genes <- read.table(file.path(chipster.tools.path, 'CanGEM', 'Ensembl_Genes', paste(genome.build, '.txt', sep='')), sep='\t', header=TRUE, row.names=1, as.is=TRUE)
colnames(genes) <- tolower(colnames(genes))
colnames(genes)[colnames(genes)=='chr'] <- 'chromosome'
colnames(genes)[colnames(genes)=='band'] <- 'cytoband'
colnames(genes)[colnames(genes)=='cnv.per.mb'] <- 'cnv.per.Mb'

# remove genes from chromosomes not present in the array data
genes <- genes[genes$chromosome %in% dat$chromosome,]

# define the functions for calculating gene copy numbers from probe copy numbers
unambiguous <- function(values) {
  values <- values[!is.na(values)]
  uniques <- unique(values)
  if (length(uniques) == 1)
    return(uniques)
  return(0)
}
majority <- function(values) {
  values <- values[!is.na(values)]
  amps <- mean(values==2)
  if (amps > 0.5) return(2)
  # treat amplifications as gains
  values[values==2] <- 1
  gains <- mean(values==1)
  if (gains > 0.5) return(1)
  losses <- mean(values==-1)
  if (losses > 0.5) return(-1)
  return(0)
}

calls <- dat[,grep("^flag\\.", colnames(dat))]
logratios <- dat[,grep("^chip\\.", colnames(dat))]

# identify additonal matrices (segmented, ...) present in the data
x <- colnames(dat)
suffix <- sub('^chip\\.', '', x[grep('^chip\\.', x)[1]])
matrices <- sub(suffix, '', x[grep(suffix, x)])
matrices <- matrices[matrices != 'flag.']
matrices <- matrices[matrices != 'chip.']
for (m in matrices)
  logratios <- cbind(logratios, dat[,grep(m, colnames(dat))])

get.gene.data <- function(x) {
  chr <- x['chromosome']
  start <- as.integer(x['start'])
  end <- as.integer(x['end'])
  # are there probes overlapping with the position of the gene
  overlapping.probes <- which(dat$chromosome == chr &
                                     dat$end >= start &
                                   dat$start <= end)
  if (length(overlapping.probes) > 0) {
    # if yes, use those probes to calculate the copy number
    gene.calls <- apply(calls[overlapping.probes,], 2, method.for.calls)
    gene.logratios <- apply(logratios[overlapping.probes,], 2, method.for.others, na.rm=TRUE)
  } else {
    # if not, use the last preceding and the first tailing probe
    preceding.probes <- which(dat$chromosome == chr &
                                     dat$end <  start)
    tailing.probes <- which(dat$chromosome == chr &
                                 dat$start >  end)
    adjacent.probes <- preceding.probes[length(preceding.probes)]
    if (length(tailing.probes) > 0)
      adjacent.probes <- c(adjacent.probes, tailing.probes[1])
    if (length(adjacent.probes) > 0) {
      gene.calls <- apply(calls[adjacent.probes,], 2, method.for.calls)
      gene.logratios <- apply(logratios[adjacent.probes,], 2, method.for.others, na.rm=TRUE)
    } else {
      gene.calls <- 0
      gene.logratios <- 0
    }
  }
  c(gene.calls, gene.logratios)
}

# first try parallel computing
prob <- TRUE
try({
  library(snowfall)
  sfInit(parallel=TRUE, cpus=4)
  sfExport(list=c('dat', 'calls', 'logratios', 'method.for.calls', 'method.for.others', 'unambiguous', 'majority'))
  gene.calls.and.logratios <- t(sfApply(genes, 1, get.gene.data))
  sfStop()
  prob <- FALSE
}, silent=TRUE)
# if problems, fall back to sequential computing
if (prob)
  gene.calls.and.logratios <- t(apply(genes, 1, get.gene.data))

calls.bygene <- gene.calls.and.logratios[,1:ncol(calls)]
genes$loss.freq <- mean(as.data.frame(t(calls.bygene==-1)))
genes$gain.freq <- mean(as.data.frame(t(calls.bygene==1)))
if (2 %in% calls.bygene)
    genes$amp.freq <- mean(as.data.frame(t(calls.bygene==2)))
genes <- cbind(genes, gene.calls.and.logratios)

write.table(genes, file='gene-aberrations.tsv', sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# EOF
