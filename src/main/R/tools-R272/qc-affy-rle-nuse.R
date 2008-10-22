# ANALYSIS "Quality control"/"Affymetrix - using RLE and NUSE" (Affymetrix quality control using NUSE and RLE. 
# This tool should be run on RAW data, i.e., CEL-files.)
# INPUT AFFY microarray[...].cel OUTPUT rle-plot.png, nuse-plot.png
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Affymetrix quality control
# JTT 9.6.2006

# Loading the libraries
library(affy)
library(affyPLM)

# Renaming variables
w<-image.width
h<-image.height

# Reading in data
dat<-ReadAffy()

# Calculating quality control values
aqc<-fitPLM(dat)

# Plotting the QC-values
par(mar=c(7, 4, 4, 2) + 0.1)
bitmap(file="rle-plot.png", width=w/72, height=h/72)
Mbox(aqc, main="RLE", las=2)
dev.off()

par(mar=c(7, 4, 4, 2) + 0.1)
bitmap(file="nuse-plot.png", width=w/72, height=h/72)
boxplot(aqc, main="NUSE", las=2)
dev.off()

