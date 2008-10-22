# ANALYSIS Preprocessing/"Remove missing values" (Removal of missing values. All observations, i.e., genes that have
# at least one missing value are excluded from the data set.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT na-omitted.tsv


# Removal of missing values
# JTT 22.6.2006

# Loads the file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Removes missing values
dat<-na.omit(dat)

# Writes a table 
write.table(dat, file=("na-omitted.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
