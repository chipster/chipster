# ANALYSIS Utilities/"Search by genomic coordinates" (Search by chromosome name and starting and ending base pair positions. The data must contain the corresponding columns.)
# INPUT GENE_EXPRS normalized.tsv
# OUTPUT search-coordinates.tsv
# PARAMETER position STRING (Position to search for. Must contain three values that are separated by tabs, hyphens, colons or two dots (e.g. X:100-200 or 7:600..700).)
# PARAMETER include.partial.overlaps [yes, no] DEFAULT yes (Whether to include only features that are completely contained within the search window, or also partial overlaps.)

# search-coordinates.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-04-13

# load inputs
dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

pos <- gsub('\\t|:|-|\\.\\.|,', ';', position)
pos <- gsub('[^-0-9XYMTxymt;]', '', pos)
pos <- strsplit(pos, ';')

if (length(pos[[1]]) != 3)
  stop('CHIPSTER-NOTE: Unsupported format (', position, ') please use e.g. "X:100-200".')

chromosome <- pos[[1]][1]
start <- as.integer(pos[[1]][2])
end <- as.integer(pos[[1]][3])

if (include.partial.overlaps == 'no') {
  dat <- dat[!is.na(dat$chromosome) &
             dat$chromosome == chromosome &
             dat$start      >= start &
             dat$end        <= end,]
} else {
  dat <- dat[!is.na(dat$chromosome) &
             dat$chromosome == chromosome &
             dat$end        >= start &
             dat$start      <= end,]
}

# write output
write.table(dat, file='search-coordinates.tsv', quote=FALSE, sep='\t')

# EOF
