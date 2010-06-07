# ANALYSIS Preprocessing/"Filter NGS results" (Allows the user to filter the results on the basis of values in any numerical column.)
# INPUT GENERIC results.tsv OUTPUT filtered-NGS-results.tsv
# PARAMETER column COLUMN_SEL (Data column to filter by)
# PARAMETER cutoff DECIMAL FROM -10000000 TO 10000000 DEFAULT 1 (Cut-off for filtering)
# PARAMETER smaller.or.larger [equal-to, smaller-than, larger-than] DEFAULT smaller-than (Smaller, equal or larger than the cutoff is filtered.)

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

