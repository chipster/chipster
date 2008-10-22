# Histograms from a table (columns=chips, rows=genes)
# JTT, 12.7.2005

# Loads the normalized data
# This needs a parameter "file" that tells which data file should be loaded.
# Datafiles are produced by normalization scripts, and at the moment
# datafile can be any of "rma.txt", "mas5.txt", or "liwong.txt".
file<-c("rma.txt")
data<-read.table(file, header=T, sep="\t")

# Calculates values needed for drawing histograms
# Values for every histogram are stored in a list, one chip after another
# histogram$mids gives the midpoints of the classes
# histogram$counts gives the number of observations at each class
len<-length(data)
histogram<-c()
for(i in 1:len) {
h<-hist(data[,i], plot=F, br=20);
histogram<-append(histogram, h)
}

# Writes the calculated values into a file
sink("histogram-values.txt")
histogram
sink()
