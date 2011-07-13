# TOOL ngs-filter-results-column.R: "Filter by column value" (Allows the user to filter the results on the basis of values in any numerical column.)
# INPUT results.tsv: results.tsv TYPE GENERIC 
# OUTPUT filtered-NGS-results.tsv: filtered-NGS-results.tsv 
# PARAMETER column: column TYPE COLUMN_SEL (Data column to filter by)
# PARAMETER has.rownames: data with rownames TYPE [yes: yes, no: no] DEFAULT no (Specifies whether the data has unique identifiers as rownmaes or lacks them.)
# PARAMETER cutoff: cutoff TYPE DECIMAL FROM -10000000 TO 10000000 DEFAULT 1 (Cut-off for filtering)
# PARAMETER smaller.or.larger: filtering criteria TYPE [equal-to: equal-to, smaller-than: smaller-than, larger-than: larger-than, within: within, outside: outside] DEFAULT smaller-than (Smaller or larger than the cutoff is filtered.
# Use the "within" or "outside" options to filter symmmetrically around two cut-offs, useful for example when searching for 2-fold up- and down-regulated genes.)

# Filter NGS results on the basis of a continuous parameter in a column
# MG 26.5.2010

# Loads the data
file <- c("results.tsv")
if (has.rownames == "yes") {
	dat <- read.table('normalized.tsv', header=TRUE, sep="\t", row.names=1, check.names=FALSE)
}
if (has.rownames == "yes") {
	dat <- read.table(file, header=T, sep="\t") 
}

# Extract the data to a vector
f <- dat[,grep(column, colnames(dat))]

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
if (has.rownames == "yes") {
	write.table(dat2, "filtered-NGS-results.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}
if (has.rownames == "no") {
	write.table(dat2, "filtered-NGS-results.tsv", sep="\t", row.names=F, col.names=T, quote=F)
}
