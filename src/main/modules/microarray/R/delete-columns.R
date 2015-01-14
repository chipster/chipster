# TOOL delete-columns.R: "Delete columns" (Delete the specified column or columns from the data.)
# INPUT normalized.tsv TYPE GENE_EXPRS
# OUTPUT deleted.tsv
# PARAMETER column1: "Column 1" TYPE COLUMN_SEL DEFAULT EMPTY (Data column to filter by)
# PARAMETER column2: "Column 2" TYPE COLUMN_SEL DEFAULT EMPTY (Data column to filter by)
# PARAMETER deletion.type: "Deletion type" TYPE [to-right, between] DEFAULT to-right (Delete columns to the right from column1 or between the specified columns)

# JTT 10.6.2008, Created
# MK 20.09.2013, fix bugs if both columns have been defined and to-right type of deletion selected
# MK 02.01.2014, added posibility to delete all columns

# Sanity checks
if(column1=="EMPTY" & column2=="EMPTY") {
	stop("CHIPSTER-NOTE: You need to select at least one column!")
}
if((column1=="EMPTY" & deletion.type=="between") | (column2=="EMPTY" & deletion.type=="between")) {
	stop("CHIPSTER-NOTE: Check the settings! Can't delete all column in between two columns, when only one column has been selected.")
}

# Loads the data
file<-c("normalized.tsv")
dat<-as.data.frame(read.table(file, sep="\t", header=T, row.names=1, quote=""))

# Filtering
if(deletion.type=="to-right") {
	column.sel <- ifelse(column1!="EMPTY", column1, column2)
	column.num <- which(names(dat)==column.sel)
	if(column.num > 1) {
		dat2 <- dat[,1:((column.num)-1)]
	} else {
		dat2 <- as.data.frame(matrix(NA, ncol=1, nrow=nrow(dat)))		
		rownames(dat2) <- rownames(dat)
		dat2 <- dat2[,0:0]
	}
} else {
	column.min <- min(which(names(dat)==column1), which(names(dat)==column2))
	column.max <- max(which(names(dat)==column1), which(names(dat)==column2))
	if(column.min==1 & column.max == ncol(dat)) {
		dat2 <- as.data.frame(matrix(NA, ncol=1, nrow=nrow(dat)))		
		rownames(dat2) <- rownames(dat)
		dat2 <- dat2[,0:0]
	} else {
		dat2<-dat[,-(column.min:column.max)]
	}
}

# Writes out the combined table
if(ncol(dat2) >= 1) {
	for(i in 1:ncol(dat2)) {
		dat2[,i] <- gsub("'", "", dat2[,i])
		dat2[,i] <- gsub("\"", "", dat2[,i])
	}
}

write.table(dat2, "deleted.tsv", sep="\t", row.names=T, col.names=T, quote=F)
