# TOOL ngs-filter-results-column.R: "Filter table by column value" (Filters a data table based on values in any numerical column.)
# INPUT results.tsv: results.tsv TYPE GENERIC 
# OUTPUT filtered-NGS-results.tsv: filtered-NGS-results.tsv 
# PARAMETER column: "Column to filter by" TYPE COLUMN_SEL (Data column to filter by)
# PARAMETER has.rownames: "Does the first column have a title" TYPE [yes: no, no: yes] DEFAULT no (Specifies whether the first column has a title or not.)
# PARAMETER cutoff: Cutoff TYPE DECIMAL FROM -10000000 TO 10000000 DEFAULT 0.05 (Cut-off for filtering)
# PARAMETER smaller.or.larger: "Filtering criteria" TYPE [equal-to: equal-to, smaller-than: smaller-than, larger-than: larger-than, within: within, outside: outside] DEFAULT smaller-than (Smaller or larger than the cutoff is filtered.
# Use the "within" or "outside" options to filter symmmetrically around two cut-offs, useful for example when searching for 2-fold up- and down-regulated genes.)

# Filter NGS results on the basis of a continuous parameter in a column
# 26.05.2010, MG 
# 08.05.2013, EK clarified wording
# 08.05.2014, MK read.table command fixed and change the script to support duplicate row names (e.g. BED files)
# 01.06.2014, EK reverted the has.rownames parameter definition, moved BED support to a separate script

# Loads the data
file <- c("results.tsv")
if (has.rownames == "yes") {
	dat <- read.table(file, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
}
if (has.rownames == "no") {
	dat <- read.table(file, header=TRUE, sep="\t", row.names=NULL) 
}

# Extract the data to a vector
f <- dat[,colnames(dat) %in% column]
if(length(which(colnames(dat) %in% column == TRUE)) > 1) {
stop("CHIPSTER-NOTE: Please choose a column that occurs only once in the table")
	}

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
