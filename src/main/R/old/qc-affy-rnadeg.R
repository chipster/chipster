# RNA degradation plots, used for quality control of Affymetrix chips
# JTT, 12.7.2005

# loads the needed libraries
library(affy)

# Reads raw data
data<-ReadAffy()

# Calculation of RNA-deg values
# deg$means.by.number gives the degradation values at each sample for every QC-probe (11)
# These values can be used for plotting the degradation values outside R 
deg<-AffyRNAdeg(data, log.it=T)

# Plotting the values into a jpeg file "RNA-degradation-plot.jpg"
jpeg(file.path(cwd, "RNA-degradation-plot.jpg"), quality=100)
plotAffyRNAdeg(deg)
dev.off()

# Degradation values can be saved to file with the following three lines
# sink("RNA-degradation-values.txt")
# deg
# sink()
# This writes the output to stdout
# dput(deg, file="")
