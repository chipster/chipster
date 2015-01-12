# TOOL cna-count-aberrations.R: "Count the number of aberrations per sample" (For each sample, counts the number of aberrations that are at least as big as the minimum size provided.)
# INPUT regions.tsv: regions.tsv TYPE GENERIC
# OUTPUT aberration-counts.tsv: aberration-counts.tsv
# PARAMETER minimum.size: "Minimum size" TYPE INTEGER FROM 0 DEFAULT 0 (The minimum aberration size to count, in base pairs.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-23

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-CGHregions.R'))

input <- readData('regions.tsv')
regions <- toCghRegions(input)

fuse.regions <- function(cgh) {
  x <- fData(cgh)[, 1:3]
  if (class(cgh) == 'cghRegions') {
    x$Call=regions(cgh)[, 1]
  } else {
    x$Call=calls(cgh)[, 1]
  }
  colnames(x)[4] <- 'Call'
  fused.data <- data.frame(Chromosome=integer(), Start=integer(), End=integer(), Call=integer())
  for (chr in unique(x$Chromosome)) {
    chr.data <- x[x$Chromosome == chr, c('Start', 'End', 'Call')]
    prev.call <- chr.data[1, 'Call']
    start <- chr.data[1, 'Start']
    if (nrow(chr.data) > 1) {
      for (i in 2:nrow(chr.data)) {
        if (chr.data[i, 'Call'] != prev.call) {
          fused.data <- rbind(fused.data, data.frame(Chromosome=chr, Start=start, End=chr.data[i-1, 'End'], Call=prev.call))
          prev.call <- chr.data[i, 'Call']
          start <- chr.data[i, 'Start']
        }
      }
      fused.data <- rbind(fused.data, data.frame(Chromosome=chr, Start=start, End=chr.data[i, 'End'], Call=prev.call))
    } else {
      fused.data <- rbind(fused.data, data.frame(Chromosome=chr, Start=start, End=chr.data[1, 'End'], Call=prev.call))
    }
  }
  rownames(fused.data) <- sprintf('%s:%i-%i', fused.data$Chromosome, fused.data$Start, fused.data$End)
  fused.data
}

m <- matrix(nrow=5, ncol=3 + ncol(regions), dimnames=list(c('deletions', 'losses', 'gains', 'amplifications', 'total'), c("min", "max", "mean", sampleNames(regions))))
for (s in sampleNames(regions)) {
  abes <- fuse.regions(regions[, s])
  abes <- abes[abes$Call != 0,]
  abes <- abes[(abes$End - abes$Start + 1) >= minimum.size,]
  m['deletions', s] <- sum(abes$Call == -2)
  m['losses', s] <- sum(abes$Call == -1)
  m['gains', s] <- sum(abes$Call == 1)
  m['amplifications', s] <- sum(abes$Call == 2)
  m['total', s] <- nrow(abes)
}
if (sum(m['amplifications', ], na.rm=TRUE) == 0)
  m <- m[-4,]
if (sum(m['deletions', ], na.rm=TRUE) == 0)
  m <- m[-1,]
m[, "min"] <- apply(m[, -(1:3)], 1, min)
m[, "max"] <- apply(m[, -(1:3)], 1, max)
m[, "mean"] <- round(apply(m[, -(1:3)], 1, mean), digits=3)

writeData(t(m), file="aberration-counts.tsv")

# EOF
