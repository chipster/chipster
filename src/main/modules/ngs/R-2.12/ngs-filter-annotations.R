# TOOL ngs-filter-annotations.R: "Filter table by column term" (Allows the user to filter the table rows on the basis of terms in any text column.)
# INPUT annotations.tsv: annotations.tsv TYPE GENERIC 
# OUTPUT filtered-NGS-results.tsv: filtered-NGS-results.tsv 
# PARAMETER column: "Column to filter by" TYPE COLUMN_SEL (Data column to filter by)
# PARAMETER match.term: "Term to match" TYPE STRING DEFAULT empty (Textual term to search for.)
# PARAMETER has.rownames: "Does the first column have a title" TYPE [yes: no, no: yes] DEFAULT no (Specifies whether the data has unique identifiers as rownames or lacks them.)


# MG 29.5.2010
# MK, EK 21.08.2013 added support for rownames


# Loads the normalized data
file<-c("annotations.tsv")
dat <- read.table(file, header=T, sep="\t", comment.char="")


if(column == " ") {
	dat2 <- dat[grep(match.term, rownames(dat)),]
} else {
# Extract the data from the column in question
	dat2 <- dat[grep(match.term, dat[,(as.vector(grep(column,names(dat))))]),]
}


# Write the data to disk
if (has.rownames == "yes") {
	write.table(dat2, "filtered-NGS-results.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}
if (has.rownames == "no") {
	write.table(dat2, "filtered-NGS-results.tsv", sep="\t", row.names=F, col.names=T, quote=F)
}