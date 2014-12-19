# TOOL sort_and_select.R: "Sort table by column value" (Sort a data table based on selected column.)
# INPUT results.tsv: results.tsv TYPE GENERIC 
# OUTPUT sorted-results.tsv: sorted-results.tsv 
# PARAMETER column: "Column to sort by" TYPE COLUMN_SEL (Data column to filter by)
# PARAMETER has.rownames: "Does the first column have a title" TYPE [yes: no, no: yes] DEFAULT no (Specifies whether the first column has a title or not.)

# Sort NGS results on the basis of a continuous parameter in a column
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


# sort the data
dat.sorted <- dat[order(column) , ]


# Writing the data to disk
if (has.rownames == "yes") {
	write.table(dat.sorted, "filtered-NGS-results.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}
if (has.rownames == "no") {
	write.table(dat.sorted, "filtered-NGS-results.tsv", sep="\t", row.names=F, col.names=T, quote=F)
}

  