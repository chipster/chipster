# TOOL search-coordinates.R: "Search by genomic coordinates" (Search by chromosome name and starting and ending base pair positions. The data must contain the corresponding columns.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT search-coordinates.tsv: search-coordinates.tsv 
# PARAMETER position: Position TYPE STRING DEFAULT empty (Position to search for. Must contain three values that are separated by tabs, hyphens, colons or two dots (e.g. X:100-200 or 7:600..700\).)
# PARAMETER include.partial.overlaps: "Include partial overlaps" TYPE [yes: yes, no: no] DEFAULT yes (Whether to include only features that are completely contained within the search window, or also partial overlaps.)

# Ilari Scheinin <firstname.lastname@gmail.com>: 2012-12-13
# MK, 25-10-2013, modified to accept also chr-columns. Modified so that rows without chr-information are removed


# load inputs
file <- 'normalized.tsv'
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, check.names=FALSE)

if(length(grep("^chromosome$", colnames(dat)))==1) {
	chr_string <- "chromosome"
} else if(length(grep("^chr$", colnames(dat)))==1) {
	chr_string <- "chr"
} else {
	chr_string <- "chromosome"
}

pos <- c(chr_string,'start','end')
if (length(setdiff(pos, colnames(dat)))!=0) {
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome (or chr), start, end.')
}

# parse coordinates
position <- gsub('(\t|:|-|\\.\\.)', ';', position)
items <- strsplit(position, ';')[[1]]
chromosome <- items[1]
start <- as.integer(items[2])
end <- as.integer(items[3])

# discard data from other chromosomes
chr_listing <- dat[, chr_string]
dat <- dat[chr_listing == chromosome,]
dat <- dat[-which(is.na(dat[,1])),]

if(length(which(dat$start >= start & dat$end <= end)) == 0) {
  stop('CHIPSTER-NOTE: Given genomic region yields no matches.')
}

if (include.partial.overlaps == 'no') {
  dat <- dat[dat$start >= start & dat$end <= end,]
} else {
  dat <- dat[dat$end > start & dat$start < end,]
}

# write output
options(scipen=10)
write.table(dat, file='search-coordinates.tsv', quote=FALSE, sep='\t')

# EOF
