# ANALYSIS Preprocessing/"Filter using a column term" (Allows the user to filter the genes on the basis of one column with textual terms, such as gene symbol, gene description or GO term.)
# INPUT GENERIC normalized.tsv OUTPUT column-term-filter.tsv
# PARAMETER column COLUMN_SEL (Data column to filter by)
# PARAMETER match.term STRING DEFAULT empty (String to search for)
# PARAMETER exact.match [yes,no] DEFAULT yes (Should only exact matches be retained?)

# Tool that allows the data to be filtered according to text terms in a column
# MG 9.3.2010
# IS 1.11.2010, small fix to column matching
# JTT 30.3.2010, added a parameter to control matching

# Loads the normalized data
dat <- read.table('normalized.tsv', header=TRUE, sep='\t', quote='')

# Only exact matches?
if(exact.match=="yes") {
   match.term<-paste("^", match.term, "$", sep="")
}

# grep the desired rows
dat2 <- dat[grep(match.term, dat[,column]),]

# write output
write.table(dat2, 'column-filter3.tsv', sep='\t', quote=FALSE)
# EOF


