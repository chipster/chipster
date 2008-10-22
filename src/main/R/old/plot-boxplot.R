# Boxplot of a table (column=chips, rows=genes)
# JTT, 12.7.2005

# Loads the normalized data
# This needs a parameter "file" that tells which data file should be loaded.
# Datafiles are produced by normalization scripts, and at the moment
# datafile can be any of "rma.txt", "mas5.txt", or "liwong.txt".
file<-c("rma.txt")
data<-read.table(file, header=T, sep="\t")

# Draws a boxplot

# Drawing needs a parameter "groups" that specifies the grouping of the samples
# This can be skipped if needed, "groups" is only used for coloring the boxplot
# If skipped, remove parameter "col" from from boxplot call
groups<-c(1,1,2,2,3,3)
# Vector can also be read from a file
# groups<-scan("groups", sep=",")

# This draws a boxplot into a jpeg file "boxplot.jpg"
cwd<-getwd()
jpeg(file.path(cwd, "boxplot.jpg"), quality=100)
boxplot(data, col=groups+1)
dev.off()

# This does not draw a histogram, but saves the values into a variable boxplotvalues
# boxplotvalues$stats specifies the five numbers giving the locations of horizontal lines in the boxplot
# boxplotvalues$out and boxplotvalues$group give the outliers and their groups
# boxplotvalues$names gives the names of the individual boxes or here, the names of the samples (chips)
boxplotvalues<-boxplot(data, col=groups+1, plot=F)

# If needed "boxplotvalues" can be written as file with any of the following commands
# 1. example
# dput(boxplotvalues, "boxplotvalues.txt")
# 2. example
# library(marray)
# write.list(boxplotvalues, filename="boxplotvalues.txt")
# 3. example
# sink("boxplotvalues.txt")
# boxplotvalues
# sink()
# This writes the output to stdout
# dput(boxplotvalues, file="")