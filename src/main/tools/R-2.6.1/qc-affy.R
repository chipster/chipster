# ANALYSIS "Quality control"/"Affymetrix basic" (Affymetrix quality control for RNA degradation and
# general quality parameters, such as scaling factor. This tool should be run on RAW data, i.e., CEL-files.
# You need to have at least two samples for this tool to work.)
# INPUT AFFY microarray[...].cel OUTPUT RNA-degradation-plot.png, simpleaffy-plot.png
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Affymetrix quality control
# JTT 9.6.2006

# Loading the libraries
library(affy)
library(simpleaffy)

# Renaming variables
w<-image.width
h<-image.height

# Reading in data
dat<-ReadAffy()

# Calculating quality control values
aqc<-qc(dat)

# Plotting the QC-values
bitmap(file="simpleaffy-plot.png", width=w/72, height=h/72)
plot(aqc)
dev.off()

# Checking the RNA degradation
deg<-AffyRNAdeg(dat)

# Saving the degradation result into a file
cols<-sample(colors(),nrow(pData(dat)))
bitmap(file="RNA-degradation-plot.png", width=w/72, height=h/72)
plotAffyRNAdeg(deg, col=cols)
legend(legend=sampleNames(dat), x="topleft", lty=1, cex=0.75, col=cols)
dev.off()
