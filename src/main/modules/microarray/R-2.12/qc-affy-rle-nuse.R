# TOOL qc-affy-rle-nuse.R: "Affymetrix - using RLE and NUSE" (Affymetrix quality control using NUSE and RLE. This tool should be run on RAW data, i.e., CEL-files.)
# INPUT microarray{...}.cel: microarray{...}.cel TYPE AFFY 
# OUTPUT rle-plot.pdf: rle-plot.pdf 
# OUTPUT nuse-plot.pdf: nuse-plot.pdf 
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


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
pdf(file="rle-plot.pdf", width=w/72, height=h/72)
Mbox(aqc, main="RLE", las=2)
dev.off()

par(mar=c(7, 4, 4, 2) + 0.1)
pdf(file="nuse-plot.pdf", width=w/72, height=h/72)
boxplot(aqc, main="NUSE", las=2)
dev.off()

