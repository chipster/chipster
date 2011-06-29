# TOOL ngs-filter-results-column.R: "Filter by column value" (Allows the user to filter the results on the basis of values in any numerical column.)
# INPUT results.tsv: results.tsv TYPE GENERIC 
# OUTPUT filtered-NGS-results.tsv: filtered-NGS-results.tsv 
# PARAMETER column: column TYPE COLUMN_SEL (Data column to filter by)
# PARAMETER cutoff: cutoff TYPE DECIMAL FROM -10000000 TO 10000000 DEFAULT 1 (Cut-off for filtering)
# PARAMETER smaller.or.larger: smaller.or.larger TYPE [equal-to: equal-to, smaller-than: smaller-than, larger-than: larger-than] DEFAULT smaller-than (Smaller, equal or larger than the cutoff is filtered.)

# Filter NGS results on the basis of a continuous parameter in a column
# MG 26.5.2010

# Loads the data
file <- c("results.tsv")
dat <- read.table(file, header=T, sep="\t") 

# Extract the data to a vector
f <- dat[,grep(column, colnames(dat))]

# Filters the data
if (smaller.or.larger == "equal-to") {
	dat2 <- dat[which(f == cutoff),]
}
if (smaller.or.larger == "smaller-than") {
	dat2 <- dat[which(f <= cutoff),]
}
if (smaller.or.larger == "larger-than") {
	dat2 <- dat[which(f >= cutoff),]
}

# Writing the data to disk
write.table(dat2, "filtered-NGS-results.tsv", sep="\t", row.names=F, col.names=T, quote=F)

