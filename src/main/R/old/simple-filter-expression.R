# ANALYSIS "Expression"/"Expression value filtering" (Expression value filtering)
# INPUT GENE_EXPRS normalised.txt OUTPUT expression-filter.txt

# JTT, 14.7.2005
# modified by JTT, 23.11.2005
# modified to save only gene names by JTT, 1.12.2005

# Testing some gargabe collection here, might be needed for large datasets
# elsewhere, also.

# Load the packages
library(genefilter)

# Loads the normalized data
# This needs a parameter "file" that tells which data file should be loaded.
# Datafiles are produced by normalization scripts, and at the moment
# datafile can be any of "rma.txt", "mas5.txt", "loess.txt", "qspline.txt", or "liwong.txt".
file<-c("normalised.txt")
data<-read.table(file, header=T, sep="\t")

# Filter
# Cut-offs (up and down) for over- and underexpressed genes are needed as parameters
# These should be logarithms at the moment, but could be easily changed; just take logarithm of input before the filter
# 2-fold up
up<-c(1) 
# 2-fold down
down<-c(-1) 

# This scales the data so that the mean of every column (chip) is equal to zero
# It just makes filtering by a certain cut-off a bit easier to use
# Note that this does not change the normalized data!
scaled.data<-genescale(data)

# These loops do the actual filtering
# Note that the original expression values are retained even after filtering

up.data<-scaled.data
len<-length(data)
for(i in 1:len) {
up.data[,i]<-ifelse(scaled.data[,i]>=up, data[,1], NA);
}
filter.up<-na.omit(up.data)
# Removes the object "up.data" as it is not needed here anymore
rm(up.data) 
# Garbage collection, will free some memory
gc() 

down.data<-scaled.data
for(i in 1:len) {
down.data[,i]<-ifelse(scaled.data[,i]<=down, data[,1], NA);
}
filter.down<-na.omit(down.data)
# Removes the object "down.data" as it is not needed here anymore
rm(down.data) 
# Garbage collection, will free some memory
gc() 

# Removes the object "data" as it is not needed here anymore
rm(data) 
# Garbage collection, will free some memory
gc() 

# Formatting gene lists as a nice table
# Takes the lenght of the list, and makes a vector repeating "up" n (length) times
u<-rep("up",length(row.names(filter.up)))
# Takes the lenght of the list, and makes a vector repeating "down" n (length) times
d<-rep("down",length(row.names(filter.down)))
# Makes a new vector of lists of "up" and "down" values
ud<-c(u,d)
# Combines the two filtered lists row-wise
l<-rbind(filter.up,filter.down)
# Makes a data frame from the combined gene list and "up"/"down" values
# list<-data.frame(l,ud) # older version, saves the data along with the gene names
list<-row.names(data.frame(l,ud))

# Saving the results
write.table(list, "expression-filter.txt", sep="\t", row.names=F, col.names=c('PROBESET_ID'), quote=F)