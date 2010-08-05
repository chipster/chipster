# ANALYSIS "aCGH tools (beta testing)"/"Smooth waves from normalized aCGH data" (Smooths the wavy pattern typically seen in aCGH profiles. Note that you need a separate calibration data set, which is not measured from cancer samples.)
# INPUT GENE_EXPRS normalized_tumor.tsv, GENE_EXPRS normalized_calib.tsv
# OUTPUT smoothed.tsv

# smooth-acgh.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-08-05

library(CGHcall)
library(NoWaves)

# load tumor data and preprocess to deal with missing values
dat <- read.table('normalized_tumor.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
dat <- data.frame(probe=rownames(dat), dat, stringsAsFactors=FALSE)
dat$chromosome[dat$chromosome=='X'] <- 23
dat$chromosome[dat$chromosome=='Y'] <- 24
dat$chromosome[dat$chromosome=='MT'] <- 25
dat$chromosome <- as.integer(dat$chromosome)
cgh <- make_cghRaw(dat)
cgh <- preprocess(cgh, nchrom=23)
cgh <- data.frame(Probe=rownames(cgh@featureData@data), cgh@featureData@data, assayDataElement(cgh, 'copynumber'))

# load calibration data and preprocess to deal with missing values
dat2 <- read.table('normalized_calib.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
dat2 <- data.frame(probe=rownames(dat2), dat2, stringsAsFactors=FALSE)
dat2$chromosome[dat2$chromosome=='X'] <- 23
dat2$chromosome[dat2$chromosome=='Y'] <- 24
dat2$chromosome[dat2$chromosome=='MT'] <- 25
dat2$chromosome <- as.integer(dat2$chromosome)
calib <- make_cghRaw(calib)
calib <- preprocess(calib, nchrom=23)
calib <- data.frame(Probe=rownames(calib@featureData@data), calib@featureData@data, assayDataElement(calib, 'copynumber'))

# dewave tumor data
library(NoWaves)
calib <- SmoothNormals(calib)
dewaved <- CorrectTumors(cgh, calib)

rownames(dewaved) <- dewaved[,1]
dewaved <- dewaved[,-1]
colnames(dewaved) <- colnames(cgh)
dewaved$chromosome <- as.character(dewaved$chromosome)
dewaved$chromosome[dewaved$chromosome=='23'] <- 'X'
dewaved$chromosome[dewaved$chromosome=='24'] <- 'Y'
dewaved$chromosome[dewaved$chromosome=='25'] <- 'MT'

dewaved[,-(1:4)] <- round(dewaved[,-(1:4)], digits=2)

write.table(dewaved, file='smoothed.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF