# TOOL qc-affy-exon-rle-nuse.R: "Affymetrix exon arrays - using RLE and NUSE" (Affymetrix quality control using NUSE and RLE. This tool should be run on RAW data, i.e., CEL-files, for exon arrays. The chip type and the summary feature have to be specified.)
# INPUT microarray{...}.cel: microarray{...}.cel TYPE AFFY 
# OUTPUT rle-plot.png: rle-plot.png 
# OUTPUT nuse-plot.png: nuse-plot.png 
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)
# PARAMETER chiptype: chiptype TYPE [empty: empty, human: human, mouse: mouse, rat: rat] DEFAULT empty ()
# PARAMETER summary.feature: summary.feature TYPE [gene: gene, exon: exon] DEFAULT gene (Output summary type)

# Affymetrix quality control
# MG 12.1.2010

# Loading the libraries
library(affy)
library(affyPLM)

# Renaming variables
w<-image.width
h<-image.height

# Reading in data
dat<-ReadAffy()

# Set up proper cdf package information

if(chiptype=="empty") {
	stop("You need to specify the chiptype. Please run the script again.")
}
if(chiptype=="human" & summary.feature=="exon") {
	dat@cdfName<-"exon.pmcdf"
	dat@annotation<-"exon.pmcdf"
}
if(chiptype=="mouse" & summary.feature=="exon") {
	dat@cdfName<-"mouseexonpmcdf"
	dat@annotation<-"mouseexonpmcdf"
}
if(chiptype=="rat" & summary.feature=="exon") {
	dat@cdfName<-"ratexonpmcdf"
	dat@annotation<-"ratexonpmcdf"
}

if(chiptype=="human" & summary.feature=="gene") {
	dat@cdfName<-"huex10stv2hsentrezgcdf"
	dat@annotation<-"huex10stv2hsentrezg.db"
}
if(chiptype=="mouse" & summary.feature=="gene") {
	dat@cdfName<-"moex10stv1mmentrezgcdf"
	dat@annotation<-"moex10stv1mmentrezg.db"
}
if(chiptype=="rat" & summary.feature=="gene") {
	dat@cdfName<-"raex10stv1rnentrezgcdf"
	dat@annotation<-"raex10stv1rnentrezg.db"
}
chiptype<-dat@annotation


# Calculating quality control values
aqc<-fitPLM(dat)

# Plotting the QC-values
par(mar=c(7, 4, 4, 2) + 0.1)
title <- paste("RLE\n(", summary.feature, " level)", sep="")
bitmap(file="rle-plot.png", width=w/72, height=h/72)
Mbox(aqc, main=title, las=2)
dev.off()

par(mar=c(7, 4, 4, 2) + 0.1)
title <- paste("NUSE\n(", summary.feature, " level)", sep="")
bitmap(file="nuse-plot.png", width=w/72, height=h/72)
boxplot(aqc, main=title, las=2)
dev.off()

