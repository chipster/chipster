# TOOL filter-by-column-value.R: "Filter using a column value" (Filters genes based on one numerical column.)
# INPUT normalized.tsv: normalized.tsv TYPE GENERIC 
# OUTPUT column-value-filter.tsv: column-value-filter.tsv
# PARAMETER column: "Column" TYPE COLUMN_SEL (Data column to filter by)
# PARAMETER cutoff: "Cut-off" TYPE DECIMAL DEFAULT 1 (Cut-off for filtering)
# PARAMETER smaller.or.larger: "Smaller or larger" TYPE [equal-to: equal-to, smaller-than: smaller-than, larger-than: larger-than, within: within, outside: outside] DEFAULT smaller-than (Smaller or larger than the cutoff is filtered.
# Use the "within" or "outside" options to filter symmmetrically around two cut-offs, useful for example when searching for 2-fold up- and down-regulated genes.)

# Filter on the basis of continuous values in a column
# JTT 31.1.2008
#
# MG, 2.3.2010 added the option to filter "withn" and "outside" a symmetrical
# range of values
# IS, 12.10.2010 removed the restrictive range from the cutoff parameter
# IS, 7.4.2011 fixed a bug where certain column names (e.g. 'p-value' did not work)
# modified, IS, 12.10.2012, to cope with tables with gene descriptions (that typically contain 's)

# Loads the normalized data
file <- 'normalized.tsv'
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, check.names=FALSE)

# Extract the data to a vector
f <- dat[, column]
# Filters the data
if(smaller.or.larger=="equal-to") {
	dat2<-dat[which(f==cutoff),]
}
if(smaller.or.larger=="smaller-than") {
	dat2<-dat[which(f<=cutoff),]
}
if(smaller.or.larger=="larger-than") {
	dat2<-dat[which(f>=cutoff),]
}
if(smaller.or.larger=="outside") {
	cutoff_2 <- -cutoff
	dat2<-dat[which(f>=cutoff | f<=cutoff_2),]
}
if(smaller.or.larger=="within") {
	cutoff_2 <- -cutoff
	dat2<-dat[-which(f>=cutoff | f <=cutoff_2),]
}

# Writing the data to disk
options(scipen=10)
write.table(dat2, "column-value-filter.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# EOF
