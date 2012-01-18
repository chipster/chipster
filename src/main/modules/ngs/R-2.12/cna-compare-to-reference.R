# TOOL cna-compare-to-reference.R: "Compare to reference sample" (This tool compares binned read counts of a test sample to binned read counts of a reference sample by taking a ratio between the two.)
# INPUT binned-hits-test.tsv: "Binned read counts of test sample" TYPE GENERIC
# INPUT binned-hits-reference.tsv: "Binned read counts of reference sample" TYPE GENERIC
# OUTPUT binned-hits-ratio.tsv: "Ratio of binned read counts"

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-12-22

test <- read.table('binned-hits-test.tsv', header=TRUE, sep='\t', row.names=1, colClasses=c('character', 'character', 'integer', 'integer', 'numeric', 'numeric'))
reference <- read.table('binned-hits-reference.tsv', header=TRUE, sep='\t', row.names=1, colClasses=c('character', 'NULL', 'NULL', 'NULL', 'numeric', 'numeric'))

# to prevent negative values
test[,4] <- test[,4] - min(test[,4], na.rm=TRUE) + 1
test[,5] <- test[,5] - min(test[,5], na.rm=TRUE) + 1
reference[,1] <- reference[,1] - min(reference[,1], na.rm=TRUE) + 1
reference[,2] <- reference[,2] - min(reference[,2], na.rm=TRUE) + 1

if (nrow(test) != nrow(reference))
  stop("CHIPSTER-NOTE: Both input files need to be binned using the same bin size.")

ratio <- test
ratio[,4] <- round(test[,4] / reference[,1], digits=2)
ratio[,5] <- round(test[,5] / reference[,2], digits=2)

ratio$count[is.nan(ratio$count)] <- NA
ratio$count[!is.na(ratio$count) & ratio$count ==-Inf] <--10
ratio$count[!is.na(ratio$count) & ratio$count == Inf] <- 10

options(scipen=10)
write.table(ratio, 'binned-hits-ratio.tsv', quote=FALSE, sep='\t', na='')

# EOF
