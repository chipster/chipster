# TOOL acgh-smooth.R: "Smooth waves from normalized aCGH data" (Smooths the wavy pattern typically seen in aCGH profiles. Note that you need a separate calibration data set, which is not measured from cancer samples. It should be measured with the same array platform and preprocessed with similar normalization settings.)
# INPUT normalized_tumor.tsv: "cancer data set" TYPE GENE_EXPRS 
# INPUT normalized_calib.tsv: "calibration data set" TYPE GENE_EXPRS 
# OUTPUT smoothed.tsv: smoothed.tsv 

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-06-27

source(file.path(chipster.common.path, "library-Chipster.R"))
source(file.path(chipster.common.path, "library-CGHcall.R"))
library(NoWaves)

# load tumor data and preprocess to deal with missing values
input1 <- readData("normalized_tumor.tsv")
cgh <- toCgh(input1, level="copynumber")
cgh <- preprocess(cgh, nchrom=23)
cgh <- data.frame(Probe=rownames(fData(cgh)), fData(cgh), copynumber(cgh), stringsAsFactors=FALSE)

# load calibration data and preprocess to deal with missing values
input2 <- readData("normalized_calib.tsv")
calib <- toCgh(input2, level="copynumber")
calib <- preprocess(calib, nchrom=23)
calib <- data.frame(Probe=rownames(fData(calib)), fData(calib), copynumber(calib), stringsAsFactors=FALSE)

# dewave tumor data
calib <- SmoothNormals(calib)
dewaved <- CorrectTumors(cgh, calib)

# format for output
rownames(dewaved) <- dewaved[,1]
dewaved <- dewaved[,-1]
colnames(dewaved)[1:3] <- c("chromosome", "start", "end")
colnames(dewaved)[-(1:3)] <- paste0("chip.", colnames(dewaved)[-(1:3)])
dewaved$chromosome <- chromosomeToCharacter(dewaved$chromosome)
dewaved[,-(1:3)] <- signif(dewaved[,-(1:3)], digits=3)
output <- addAnnotationColumns(input1, dewaved)

writeData(output, "smoothed.tsv")

# EOF
