# ANALYSIS "aCGH tools (beta testing)"/"Convert called aCGH data from probes to genes" (Using the chromosomal locations of the probes of an aCGH data set, convert the data from probe-based to gene-based. The probes from which the copy number call for a given gene are determined as follows. If there are probes overlapping with the position of the gene, they are used. In case of no overlaps, the last preceding and first tailing probe are used.)
# INPUT GENE_EXPRS aberrations.tsv
# OUTPUT gene-aberrations.tsv
# PARAMETER call.method [majority, unambiguous] DEFAULT majority (The method majority means that if more than 50% of these probes give an aberrated signal, that call is used for the gene. The unambiguous method requires that all of the probes have the same call, otherwise the gene will be labeled as normal.)
# PARAMETER logratio.method [mean, median] DEFAULT mean (The probes from which the copy number call for a given gene are determined as follows. If there are probes overlapping with the position of the gene, they are used. In case of no overlaps, the last preceding and first tailing probe are used. This parameter specifies whether to use the mean or the median for calculating log ratios for the genes.)
# PARAMETER genome.build [GRCh37, NCBI36, NCBI35, NCBI34] DEFAULT GRCh37 (The genome build to use for fetching the gene coordinates.)

# convert-cn-probes-to-genes.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-08-16

dat <- read.table('aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

# load genes
genes <- read.table(paste('http://www.cangem.org/download.php?platform=CG-PLM-26&flag=', genome.build, sep=''), sep='\t', header=TRUE, row.names=1, as.is=TRUE)
colnames(genes) <- tolower(colnames(genes))
colnames(genes)[colnames(genes)=='chr'] <- 'chromosome'

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

calls <- dat[,grep("flag", names(dat))]
logratios <- dat[,grep("chip", names(dat))]

calls.bygene <- matrix(nrow=nrow(genes), ncol=ncol(calls), dimnames=list(rownames(genes), colnames(calls)))
logratios.bygene <- matrix(nrow=nrow(genes), ncol=ncol(logratios), dimnames=list(rownames(genes), colnames(logratios)))

for (gene in rownames(genes)) {
  # are there probes overlapping with the position of the gene
  overlapping.probes <- which(dat$chromosome == genes[gene, 'chromosome'] &
                                     dat$end >= genes[gene, 'start'] &
                                   dat$start <= genes[gene, 'end'])
  if (length(overlapping.probes) > 0) {
    # if yes, use those probes to calculate the copy number
    calls.bygene[gene,] <- apply(calls[overlapping.probes,], 2, call.method)
    logratios.bygene[gene,] <- apply(logratios[overlapping.probes,], 2, logratio.method, na.rm=TRUE)
  } else {
    # if not, use the last preceding and the first tailing probe
    preceding.probes <- which(dat$chromosome == genes[gene, 'chromosome'] &
                                     dat$end <  genes[gene, 'start'])
    tailing.probes <- which(dat$chromosome == genes[gene, 'chromosome'] &
                                 dat$start >  genes[gene, 'end'])
    adjacent.probes <- preceding.probes[length(preceding.probes)]
    if (length(tailing.probes) > 0)
      adjacent.probes <- c(adjacent.probes, tailing.probes[1])
    if (length(adjacent.probes) > 0) {
      calls.bygene[gene,] <- apply(calls[adjacent.probes,], 2, call.method)
      logratios.bygene[gene,] <- apply(logratios[adjacent.probes,], 2, logratio.method, na.rm=TRUE)
    } else {
      calls.bygene[gene,] <- 0
      logratios.bygene[gene,] <- 0
    }
  }
}

genes$loss.freq <- mean(as.data.frame(t(calls.bygene==-1)))
genes$gain.freq <- mean(as.data.frame(t(calls.bygene==1)))
if (2 %in% calls.bygene)
    genes$amp.freq <- mean(as.data.frame(t(calls.bygene==2)))
genes <- cbind(genes, calls.bygene)
genes <- cbind(genes, logratios.bygene)

write.table(genes, file='gene-aberrations.tsv', sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# EOF