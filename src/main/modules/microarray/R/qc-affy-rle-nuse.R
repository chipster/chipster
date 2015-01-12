# TOOL qc-affy-rle-nuse.R: "Affymetrix - using RLE and NUSE" (Affymetrix quality control using NUSE and RLE. This tool should be run on RAW data, i.e., CEL-files.)
# INPUT microarray{...}.cel: microarray{...}.cel TYPE AFFY 
# OUTPUT rle-plot.pdf: rle-plot.pdf 
# OUTPUT nuse-plot.pdf: nuse-plot.pdf 
# PARAMETER custom.chiptype: "Custom chiptype" TYPE [empty: empty, hgu133ahsentrezg(hgu133a): hgu133ahsentrezg(hgu133a), hgu133a2hsentrezg(hgu133av2): hgu133a2hsentrezg(hgu133av2), hgu133phsentrezg(hgu133plus): hgu133phsentrezg(hgu133plus), hgu133plus2hsentrezg(hgu133plus2): hgu133plus2hsentrezg(hgu133plus2), hgu133bhsentrezg(hgu133b): hgu133bhsentrezg(hgu133b), hgu95av2hsentrezg(hgu95av2): hgu95av2hsentrezg(hgu95av2), moe430ammentrezg(moe430a): moe430ammentrezg(moe430a), moe430bmmentrezg(moe430b): moe430bmmentrezg(moe430b), mouse430a2mmentrezg(mouse430a2): mouse430a2mmentrezg(mouse430a2), mouse4302mmentrezg(mouse4302): mouse4302mmentrezg(mouse4302), mm74av1mmentrezg(mgu74a): mm74av1mmentrezg(mgu74a), mgu74av2mmentrezg(mgu74av2): mgu74av2mmentrezg(mgu74av2), mgu74bv2mmentrezg(mgu74bv2): mgu74bv2mmentrezg(mgu74bv2), mgu74cv2mmentrezg(mgu74cv2): mgu74cv2mmentrezg(mgu74cv2), rae230arnentrezg(rae230a): rae230arnentrezg(rae230a), rae230brnentrezg(rae230b): rae230brnentrezg(rae230b), rat2302rnentrezg(rat2302): rat2302rnentrezg(rat2302), rgu34arnentrezg(rgu34a): rgu34arnentrezg(rgu34a), rgu34brnentrezg(rgu34b): rgu34brnentrezg(rgu34b), rgu34crnentrezg(rgu34c): rgu34crnentrezg(rgu34c)] DEFAULT empty (Custom chiptype. If not given, inferred from the chips themselves)
# PARAMETER image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# Affymetrix quality control
# JTT 9.6.2006
# MK: 12.06.2013 added possibility to analyse custom chips

# Loading the libraries
library(affy)
library(affyPLM)

# Renaming variables
w<-image.width
h<-image.height

# Reading in data
dat<-ReadAffy()

# Use custom annotations
if(custom.chiptype!="empty") {
	chiptype <- custom.chiptype
	chiptype <- gsub("\\(.*?\\)", "", chiptype) 
	dat@annotation<-chiptype
	dat@cdfName<-chiptype
}

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

