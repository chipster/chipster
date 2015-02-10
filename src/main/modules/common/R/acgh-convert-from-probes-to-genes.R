# TOOL acgh-convert-from-probes-to-genes.R: "Detect genes from called copy number data" (Using the chromosomal locations of the data points of a copy number data set, convert the data from probe/bin-based to gene-based. The probes/bins from which the copy number call for a given gene are determined as follows. If there are probes/bins overlapping with the position of the gene, they are used. In case of no overlaps, the last preceding and first tailing probe/bin are used.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# OUTPUT gene-aberrations.tsv: gene-aberrations.tsv 
# PARAMETER method.for.calls: "Method for calls" TYPE [majority: majority, unambiguous: unambiguous] DEFAULT majority (The method majority means that if more than 50% of these probes/bins give an aberrated signal, that call is used for the gene. The unambiguous method requires that all of the probes/bins have the same call, otherwise the gene will be labeled as normal.)
# PARAMETER method.for.others: "Method for others" TYPE [mean: mean, median: median] DEFAULT mean (Whether to use the mean or the median for calculating other data than copy number calls.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-24

source(file.path(chipster.common.path, 'library-Chipster.R'))
library(Homo.sapiens)

dat <- readData("aberrations.tsv")

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

# load genes
all.genes <- AnnotationDbi::select(Homo.sapiens::Homo.sapiens, keys=keys(Homo.sapiens::Homo.sapiens,
    keytype='CHRLOC'), columns=c('CHRLOC', 'CHRLOCEND', 'ENTREZID', 'SYMBOL', 'GENENAME'), keytype='CHRLOC')
all.genes <- all.genes[all.genes$CHRLOCCHR %in% unique(dat$chromosome), ]
all.genes$CHRLOC <- abs(all.genes$CHRLOC)
all.genes$CHRLOCEND <- abs(all.genes$CHRLOCEND)
all.genes <- all.genes[, c('CHRLOCCHR', 'CHRLOC', 'CHRLOCEND', 'SYMBOL', 'GENENAME', 'ENTREZID')]
colnames(all.genes) <- c('chromosome', 'start', 'end', 'symbol', 'description', 'entrez')

all.genes$entrez.chr <- paste(all.genes$entrez, all.genes$chromosome, sep=";")
ambiguous <- unique(all.genes$entrez.chr[duplicated(all.genes$entrez.chr)])
unambiguous.positions <- all.genes[!all.genes$entrez.chr %in% ambiguous, ]
positions <- vector("list", length(ambiguous))
for (i in seq_along(ambiguous)) {
  entries <- all.genes[all.genes$entrez.chr == ambiguous[i], ]
  entries$start <- min(entries$start)
  entries$end <- max(entries$end)
  # entries <- unique(entries)
  # if (nrow(entries) > 1) stop("PERKELE!")
    positions[[i]] <- entries[1, ]
}
genes <- do.call(rbind, positions)
genes[(nrow(genes)+1):(nrow(genes)+nrow(unambiguous.positions)),] <- unambiguous.positions
genes$entrez.chr <- NULL
genes <- genes[order(chromosomeToInteger(genes$chromosome), genes$start, genes$end), ]
rownames(genes) <- 1:nrow(genes)

all.genes <- rm(list=c("all.genes", "ambiguous", "unambiguous.positions", "positions", "entries", "i"))

# define the functions for calculating gene copy numbers from probe/bin copy numbers
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
  # are there probes/bins overlapping with the position of the gene
  overlapping.probes <- which(dat$chromosome == chr &
                                     dat$end >= start &
                                   dat$start <= end)
  if (length(overlapping.probes) > 0) {
    # if yes, use those probes/bins to calculate the copy number
    gene.calls <- apply(calls[overlapping.probes,], 2, method.for.calls)
    gene.logratios <- apply(logratios[overlapping.probes,], 2, method.for.others, na.rm=TRUE)
  } else {
    # if not, use the last preceding and the first tailing probe/bin
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

calls.bygene <- gene.calls.and.logratios[, 1:ncol(calls)]
if (-2 %in% calls.bygene)
  genes$del.freq <- round(rowMeans(calls.bygene == -2), digits=3)
genes$loss.freq <- round(rowMeans(calls.bygene == -1), digits=3)
genes$gain.freq <- round(rowMeans(calls.bygene == 1), digits=3)
if (2 %in% calls.bygene)
  genes$amp.freq <- round(rowMeans(calls.bygene == 2), digits=3)
genes <- cbind(genes, gene.calls.and.logratios)

writeData(genes, "gene-aberrations.tsv")

# EOF
