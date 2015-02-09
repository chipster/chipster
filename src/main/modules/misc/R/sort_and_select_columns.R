# TOOL sort_and_select_columns.R: "Sort table by column value" (Sort a data table based on selected column.)
# INPUT results.tsv: results.tsv TYPE GENERIC 
# OUTPUT sorted-results.tsv: sorted-results.tsv 
# OUTPUT OPTIONAL sort.log
# PARAMETER column: "Column to sort by" TYPE COLUMN_SEL (Data column to sort by)
# PARAMETER sortorder: "Sorting order" TYPE [ascending: "Ascending", descending: "Descending"] DEFAULT ascending (Sortitng order)


# Load the data
file <- c("results.tsv")

dat <- read.table(file, header=TRUE, sep="\t", check.names=FALSE, quote="", stringsAsFactors=FALSE)

# sort the data
if (sortorder == "descending") {
	dat.sorted <- dat[order(dat[[column]], decreasing=TRUE), ]	
}

if (sortorder == "ascending") {
	dat.sorted <- dat[order(dat[[column]]), ]
}

write.table(dat.sorted, "sorted-results.tsv", sep="\t", row.names=F, col.names=T, quote=F)


  