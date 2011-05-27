# TOOL delete-columns.R: "Delete columns" (Delete the specified column or columns from the data.)
# INPUT normalized.tsv TYPE GENE_EXPRS
# OUTPUT deleted.tsv
# PARAMETER column1 TYPE COLUMN_SEL DEFAULT EMPTY (Data column to filter by)
# PARAMETER column2 TYPE COLUMN_SEL DEFAULT EMPTY (Data column to filter by)
# PARAMETER deletion.type TYPE [to-right, between] DEFAULT to-right (Delete columns to the right from column1 or between the specified columns)


# Deletes columns from data
# JTT 10.6.2008

# Parameter settings (default) for testing purposes
#column1<-c("chip.microarray017.cel")
#column2<-c("EMPTY")
#deletion.type<-c("to-right")

# Sanity checks
if(column1=="EMPTY" & column2=="EMPTY") {
	stop("You need to select at least one column!")
}
if((column1=="EMPTY" & deletion.type=="between") | (column2=="EMPTY" & deletion.type=="between")) {
	stop("Check the settings! Can't delete all column in between two columns, when only one column has been selected.")
}

# Loads the data
file<-c("normalized.tsv")
dat<-read.table(file, sep="\t", header=T, row.names=1)

# Filtering
if(column1!="EMPTY" & deletion.type=="to-right") {
	dat2<-dat[,1:(which(names(dat)==column1)-1)]
}

if(column2!="EMPTY" & deletion.type=="to-right") {
	dat2<-dat[,1:(which(names(dat)==column2)-1)]
}

if(column1!="EMPTY" & column2!="EMPTY" & deletion.type=="between") {
	dat2<-dat[,-(min(which(names(dat)==column1), which(names(dat)==column2)):max(which(names(dat)==column1), which(names(dat)==column2)))]
}

# Writes out the combined table
write.table(dat2, "deleted.tsv", sep="\t", row.names=T, col.names=T, quote=F)
