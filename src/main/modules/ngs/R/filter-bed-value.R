# TOOL filter-bed-value.R: "Filter BED by column value" (Filters a BED file based on a numerical column.)
# INPUT regions.bed: regions.bed TYPE GENERIC 
# OUTPUT filtered.bed: filtered.bed 
# PARAMETER column: "Column to filter by" TYPE COLUMN_SEL (Data column to filter by)
# PARAMETER cutoff: Cutoff TYPE DECIMAL FROM -10000000 TO 10000000 DEFAULT 0.05 (Cut-off for filtering)
# PARAMETER smaller.or.larger: "Filtering criteria" TYPE [equal-to: equal-to, smaller-than: smaller-than, larger-than: larger-than, within: within, outside: outside] DEFAULT smaller-than (Smaller or larger than the cutoff is filtered.
# Use the "within" or "outside" options to filter symmmetrically around two cut-offs, useful for example when searching for up- and down-regulated genes.)

# 01.06.2014, EK 

# Read the file
file <- c("regions.bed")
dat <- read.table(file, header=FALSE, sep="\t", row.names=NULL) 

# Extract the data to a vector
column <- as.numeric(gsub("column", "", column))+1
f <- dat[, column]	

# Filters the data
if(smaller.or.larger=="equal-to") {
	dat2<-dat[which(f==cutoff),]
}
if(smaller.or.larger=="smaller-than") {
	dat2<-dat[which(f<=cutoff),]
}
if(smaller.or.larger=="larger-than") {
	dat2<-dat[which(f>=cutoff),]
}
if(smaller.or.larger=="outside") {
	cutoff_2 <- -cutoff
	dat2<-dat[which(f>=cutoff | f<=cutoff_2),]
}
if(smaller.or.larger=="within") {
	cutoff_2 <- -cutoff
	dat2<-dat[-which(f>=cutoff | f <=cutoff_2),]
}

# Writing the data to disk
write.table(dat2, "filtered.bed", sep="\t", row.names=F, col.names=F, quote=F)
