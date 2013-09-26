# TOOL delete-columns.R: "Delete columns" (Delete the specified column or columns from the data.)
# INPUT normalized.tsv TYPE GENE_EXPRS
# OUTPUT deleted.tsv
# PARAMETER column1 TYPE COLUMN_SEL DEFAULT EMPTY (Data column to filter by)
# PARAMETER column2 TYPE COLUMN_SEL DEFAULT EMPTY (Data column to filter by)
# PARAMETER deletion.type TYPE [to-right, between] DEFAULT to-right (Delete columns to the right from column1 or between the specified columns)

# JTT 10.6.2008, Created
# MK 20.09.2013, fix bugs if both columns have been defined and to-right type of deletion selected

# Sanity checks
if(column1=="EMPTY" & column2=="EMPTY") {
	stop("You need to select at least one column!")
}
if((column1=="EMPTY" & deletion.type=="between") | (column2=="EMPTY" & deletion.type=="between")) {
	stop("Check the settings! Can't delete all column in between two columns, when only one column has been selected.")
}

# Loads the data
file<-c("normalized.tsv")
dat<-read.table(file, sep="\t", header=T, row.names=1, quote="")

# Filtering
if(deletion.type=="to-right") {
	column.sel <- ifelse(column1!="EMPTY", column1, column2)
	dat2<-dat[,1:(which(names(dat)==column.sel)-1)]	 
} else {
	dat2<-dat[,-(min(which(names(dat)==column1), which(names(dat)==column2)):max(which(names(dat)==column1), which(names(dat)==column2)))]
}

# Writes out the combined table
for(i in 1:ncol(dat2)) {
	dat2[,i] <- gsub("'", "", dat2[,i])
	dat2[,i] <- gsub("\"", "", dat2[,i])
}

write.table(dat2, "deleted.tsv", sep="\t", row.names=T, col.names=T, quote=F)
