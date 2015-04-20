# TOOL merge-tables.R: "Merge tables" (Merges two tables using identifiers in the first column.)
# INPUT table1.tsv: "Table 1" TYPE GENERIC
# INPUT table2.tsv: "Table 2" TYPE GENERIC
# OUTPUT combined.tsv: "Merged table"
# PARAMETER OPTIONAL include.everything: "Include all the rows in the result file" TYPE [yes, no] DEFAULT no (Include also the non-matching lines in the result file.)

# JTT 22.10.2007
# EK 20.4.2015 clarified the script

# Loads the tables
table1<-read.table(file="table1.tsv", sep="\t", header=T, row.names=1)
table2<-read.table(file="table2.tsv", sep="\t", header=T, row.names=1)

for (i in 1:ncol(table1)) {
	table1[,i] <- gsub("\t+", " ", table1[,i], perl=T)
	table1[,i] <- gsub("\n+", " ", table1[,i], perl=T)
	table1[,i] <- gsub(" +", " ", table1[,i], perl=T)
	
	table1[,i] <- gsub("\"+", ",", table1[,i], perl=T)
	table1[,i] <- gsub("\'+", ",", table1[,i], perl=T)
}

# Combines tables using row names
include=FALSE
if( include.everything == "yes" ) include=TRUE
table3<-merge(table1, table2, by.x="row.names", by.y="row.names", all=include)
row.names(table3)<-table3$Row.names
table3<-table3[-1]

# Writes out the combined table
write.table(table3, "combined.tsv", sep="\t", row.names=T, col.names=T, quote=F)
