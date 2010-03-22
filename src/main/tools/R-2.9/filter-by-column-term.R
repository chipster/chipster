# ANALYSIS Preprocessing/"Filter using a column term" (Allows the user to filter the genes on the basis of one column with textual terms, such as gene symbol, gene description or GO term.)
# INPUT GENERIC normalized.tsv OUTPUT column-filter3.tsv
# PARAMETER column COLUMN_SEL (Data column to filter by)
# PARAMETER match.term STRING DEFAULT empty (String to search for.)

# Tool that allows the data to be filtered according to text terms in a column
# MG 9.3.2010

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)
dat2<-dat[grep(match.term, dat[,(as.vector(grep(column,names(dat))))]),]
write.table(dat2, "column-filter3.tsv", sep="\t", row.names=T, col.names=T, quote=F)


