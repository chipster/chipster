# TOOL sort_and_select_columns.R: "Sort table by column" (Sort a data table based on a selected column.)
# INPUT results.tsv: results.tsv TYPE GENERIC 
# OUTPUT sorted-results.tsv: sorted-results.tsv 
# OUTPUT OPTIONAL sort.log
# PARAMETER column: "Column to sort by" TYPE COLUMN_SEL (Data column to sort by)
# PARAMETER sortorder: "Sorting order" TYPE [ascending: "Ascending", descending: "Descending"] DEFAULT ascending (Sortitng order)
# PARAMETER OPTIONAL rstyle: "First column has rownames" TYPE [yes: Yes, no: No] DEFAULT no (Choose Yes if the first column does not have a title.)

# KM 8.2.2015
# KM 17.4.2015 Added support for rownames


# Load the data
file <- c("results.tsv")

if (rstyle=="yes"){
	dat <- read.table(file, header=TRUE, sep="\t", check.names=TRUE, quote="", stringsAsFactors=FALSE)
}

if (rstyle=="no"){
	dat <- read.table(file, header=TRUE, sep="\t", check.names=FALSE, quote="", stringsAsFactors=FALSE)
}

# sort the data
if (sortorder == "descending") {
	dat.sorted <- dat[order(dat[[column]], decreasing=TRUE), ]	
}

if (sortorder == "ascending") {
	dat.sorted <- dat[order(dat[[column]]), ]
}


if (rstyle=="no"){
	write.table(dat.sorted, "sorted-results.tsv", sep="\t", row.names=F, col.names=T, quote=F)
}

if (rstyle=="yes"){
	write.table(dat.sorted, "sorted-results.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}


  