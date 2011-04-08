# ANALYSIS Preprocessing/"Filter NGS annotations" (Allows the user to filter the table rows on the basis of terms in any text column.)
# INPUT GENERIC annotations.tsv OUTPUT filtered-NGS-annotations.tsv
# PARAMETER column COLUMN_SEL (Data column to filter by)
# PARAMETER match.term STRING DEFAULT empty (String to search for.)

# Filter NGS annotation tables on the basis of a text term in a column
# MG 29.5.2010

# Loads the normalized data
file<-c("annotations.tsv")
dat <- read.table(file, header=T, sep="\t", comment.char="")

# Extract the data from the column in question
dat2 <- dat[grep(match.term, dat[,(as.vector(grep(column,names(dat))))]),]

# Write the data to disk
write.table (dat2, "filtered-NGS-annotations.tsv", sep="\t", row.names=F, col.names=T, quote=F)

