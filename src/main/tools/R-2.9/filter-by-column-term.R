# ANALYSIS Preprocessing/"Filter using a column term" (Allows the user to filter the genes on the basis of one column with textual terms, such as gene symbol, gene description or GO term.)
# INPUT GENERIC normalized.tsv OUTPUT column-filter3.tsv
# PARAMETER column COLUMN_SEL (Data column to filter by)
# PARAMETER match.term STRING DEFAULT empty (String to search for. Entering e.g. GENE1 will also match GENE10 and MYGENE1. If you want only exact matches, enter the string as ^GENE1$ - the so-called regular expression syntax.)

# Tool that allows the data to be filtered according to text terms in a column
# MG 9.3.2010
# IS 29.10.2010, small fix to column matching

# Loads the normalized data
dat <- read.table("normalized.tsv", header=TRUE, sep="\t", row.names=1)

# grep the desired rows
dat2 <- dat[grep(match.term, dat[,column]),]

# write output
write.table(dat2, "column-filter3.tsv", sep="\t", quote=FALSE)

# EOF