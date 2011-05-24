# ANALYSIS "aCGH"/"Match copy number and expression probes" (Matches the probes of a copy number data set with probes of an expression data set, using their chromosomal locations. Running this tool is a prerequisite for testing copy-number-induced effects on expression.)
# INPUT GENE_EXPRS aberrations.tsv, GENE_EXPRS normalized.tsv, GENERIC phenodata_cgh.tsv, GENERIC phenodata_exp.tsv
# OUTPUT matched-cn-and-expression.tsv, matched-cn-and-expression-heatmap.pdf, matched-phenodata.tsv
# PARAMETER sample.identifiers.1 METACOLUMN_SEL DEFAULT Sample (The phenodata column for data set 1 used to link the two data sets together.)
# PARAMETER sample.identifiers.2 METACOLUMN_SEL DEFAULT Sample (The phenodata column for data set 2 used to link the two data sets together.)
# PARAMETER method [distance, overlap, overlapplus] DEFAULT distance (The method for linking copy number and expression probes together.)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image. Not used anymore as plotting format is now PDF.)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image. Not used anymore as plotting format is now PDF.)

# match-cn-and-expression-probes.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-05-24

library(CGHcall)
library(intCNGEan)

dat1 <- read.table('aberrations.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
dat2 <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
phenodata1 <- read.table("phenodata_cgh.tsv", header=T, sep='\t', as.is=TRUE)
phenodata2 <- read.table("phenodata_exp.tsv", header=T, sep='\t', as.is=TRUE)

# determine which dataset is cgh
if (length(grep("^probnorm\\.", names(dat1)))!=0) {
  cgh <- dat1
  exp <- dat2
  phenodata_cgh <- phenodata1
  phenodata_exp <- phenodata2
  samples_cgh <- sample.identifiers.1
  samples_exp <- sample.identifiers.2
} else if (length(grep("^probnorm\\.", names(dat2)))!=0) {
  cgh <- dat2
  exp <- dat1
  phenodata_cgh <- phenodata2
  phenodata_exp <- phenodata1
  samples_cgh <- sample.identifiers.2
  samples_exp <- sample.identifiers.1
} else {
  stop('CHIPSTER-NOTE: Could not detect the aCGH data set.')
}

# check that both data sets have the probe position information
pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(cgh)))!=0)
  stop('CHIPSTER-NOTE: Both data sets must have the following columns: chromosome, start, end.')
if (length(setdiff(pos, colnames(exp)))!=0)
  stop('CHIPSTER-NOTE: Both data sets must have the following columns: chromosome, start, end.')

# check for unambiguity of sample identifiers
if (nrow(phenodata_cgh)!=length(unique(phenodata_cgh[,samples_cgh])))
  stop('CHIPSTER-NOTE: Unambigous aCGH sample identifiers: ', paste(phenodata_cgh[duplicated(phenodata_cgh[,samples_cgh]),samples_cgh], collapse=', ')) 
if (nrow(phenodata_exp)!=length(unique(phenodata_exp[,samples_exp])))
  stop('CHIPSTER-NOTE: Unambigous expression sample identifiers: ', paste(phenodata_exp[duplicated(phenodata_exp[,samples_exp]),samples_exp], collapse=', ')) 

common.samples <- intersect(phenodata_cgh[,samples_cgh], phenodata_exp[,samples_exp])
rownames(phenodata_cgh) <- phenodata_cgh[,samples_cgh]
rownames(phenodata_exp) <- phenodata_exp[,samples_exp]
phenodata_cgh$n <- 1:nrow(phenodata_cgh)
phenodata_exp$n <- 1:nrow(phenodata_exp)
index_cgh <- phenodata_cgh[common.samples, 'n']
index_exp <- phenodata_exp[common.samples, 'n']

# build a cghCall object from the cgh data
cgh$chromosome[cgh$chromosome=='X'] <- 23
cgh$chromosome[cgh$chromosome=='Y'] <- 24
cgh$chromosome[cgh$chromosome=='MT'] <- 25
cgh$chromosome <- as.integer(cgh$chromosome)
cgh <- cgh[!is.na(cgh$chromosome),]
calls <- as.matrix(cgh[,grep("^flag\\.", names(cgh))])[,index_cgh]
copynumber <- as.matrix(cgh[,grep("^chip\\.", names(cgh))])[,index_cgh]
segmented <- as.matrix(cgh[,grep("^segmented\\.", names(cgh))])[,index_cgh]
probloss <- as.matrix(cgh[,grep("^probloss\\.", names(cgh))])[,index_cgh]
probnorm <- as.matrix(cgh[,grep("^probnorm\\.", names(cgh))])[,index_cgh]
probgain <- as.matrix(cgh[,grep("^probgain\\.", names(cgh))])[,index_cgh]
if (2 %in% calls) {
  probamp <- as.matrix(cgh[,grep("^probamp\\.", names(cgh))])[,index_cgh]
  probgain <- probgain + probamp
  calls[calls==2] <- 1
}
cgh <- new('cghCall', assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=cgh$chromosome, Start=cgh$start, End=cgh$end, row.names=row.names(cgh))))
sampleNames(cgh) <- phenodata_cgh[common.samples, samples_cgh]

# build an ExpressionSet object from the expression data
dat_exp <- exp
exp$chromosome[exp$chromosome=='X'] <- 23
exp$chromosome[exp$chromosome=='Y'] <- 24
exp$chromosome[exp$chromosome=='MT'] <- 25
exp$chromosome <- as.integer(exp$chromosome)
exp <- exp[!is.na(exp$chromosome),]
exprs <- as.matrix(exp[,grep("^chip\\.", names(exp))])[,index_exp]
exp <- new("ExpressionSet", exprs=exprs, featureData=new("AnnotatedDataFrame", data=data.frame(Chromosome=exp$chromosome, Start=exp$start, End=exp$end, row.names=row.names(exp))))
sampleNames(exp) <- phenodata_exp[common.samples, samples_exp]

# match probes
matched <- intCNGEan.match(cgh, exp, CNbpend='yes', GEbpend='yes', method=method)

# bugfix version of the plotting function
intCNGEan.heatmaps.modified <- function (CNdata, GEdata, location = "mode", colorbreaks = "equiquantiles") 
{
    makeSegments <- function(data) {
        previous <- 2000
        values <- c()
        start <- c()
        end <- c()
        for (i in 1:length(data)) {
            if (data[i] != previous) {
                start <- c(start, i)
                last <- i - 1
                if (last > 0) 
                  end <- c(end, last)
                values <- c(values, data[i])
            }
            previous <- data[i]
        }
        end <- c(end, length(data))
        result <- cbind(values, start, end)
        return(result)
    }
    makeRegions <- function(CNprobs) {
        splitter <- list()
        splitter[[1]] <- c(1)
        index.temp <- 1
        j <- 1
        for (i in 1:(dim(CNprobs)[1] - 1)) {
            if (all(CNprobs[i, ] == CNprobs[i + 1, ])) {
                index.temp <- c(index.temp, i + 1)
                splitter[[j]] <- index.temp
            }
            else {
                index.temp <- i + 1
                j <- j + 1
                splitter[[j]] <- index.temp
            }
        }
        regDetails <- NULL
        for (i in 1:length(splitter)) {
            regDetails <- rbind(regDetails, c(min(splitter[[i]]), 
                max(splitter[[i]])))
        }
        return(regDetails)
    }
    if (dim(fData(CNdata))[1] != dim(fData(GEdata))[1]) {
        stop("CN and GE data have different number of rows.")
    }
    if (!all(fData(CNdata)[, 1] == fData(GEdata)[, 1])) {
        stop("chrosome annotation between CN and GE does not match.")
    }
    if (!(location %in% c("mode", "median", "mean"))) {
        stop("location parameter ill-specified.")
    }
    if (!(colorbreaks %in% c("equidistant", "equiquantiles"))) {
        stop("colorbreaks parameter ill-specified.")
    }
    exprsTemp <- as.numeric(exprs(GEdata))
    histres <- hist(as.numeric(exprs(GEdata)), plot = FALSE, 
        n = 511)
    if (location == "median") {
        exprsMode <- median(exprsTemp)
    }
    if (location == "mean") {
        exprsMode <- mean(exprsTemp)
    }
    if (location == "mode") {
        exprsMode <- histres$mids[which.max(histres$density)]
    }
    exprsTempBelowMode <- exprsTemp[exprsTemp < exprsMode]
    exprsTempAboveMode <- exprsTemp[exprsTemp >= exprsMode]
    exprsTempBelowMode <- cbind(exprsTempBelowMode, ecdf(exprsTempBelowMode)(exprsTempBelowMode))[order(exprsTempBelowMode), 
        ]
    exprsTempAboveMode <- cbind(exprsTempAboveMode, ecdf(exprsTempAboveMode)(exprsTempAboveMode))[order(exprsTempAboveMode), 
        ]
    if (colorbreaks == "equiquantiles") {
        histresBM <- hist(exprsTempBelowMode[, 2], plot = FALSE, 
            n = 100)
        histresAM <- hist(exprsTempAboveMode[, 2], plot = FALSE, 
            n = 101)
        breaks <- c(quantile(exprsTempBelowMode[, 1], probs = histresBM$breaks), 
            exprsMode, quantile(exprsTempAboveMode[, 1], probs = histresAM$breaks))
        collist <- c(maPalette(low = "red", high = "black", k = length(histresBM$breaks)), 
            maPalette(low = "black", high = "green", k = length(histresAM$breaks)))
    }
    if (colorbreaks == "equidistant") {
        collistBelowMode <- unique(maPalette(low = "red", high = "black", 
            k = 100))
        collistAboveMode <- unique(maPalette(low = "black", high = "green", 
            k = 100))
        breaks <- c(seq(min(exprsTemp), exprsMode, length.out = length(collistBelowMode) + 
            1), seq(exprsMode, max(exprsTemp), length.out = length(collistAboveMode))[-1])
        collist <- unique(c(collistBelowMode, collistAboveMode))
    }
    CNprobs <- numeric()
    for (i in 1:dim(calls(CNdata))[2]) {
        CNprobs <- cbind(CNprobs, cbind(probloss(CNdata)[, i], 
            probnorm(CNdata)[, i], probgain(CNdata)[, i], probamp(CNdata)[, 
                i]))
    }
    nclass <- dim(CNprobs)[2]/dim(calls(CNdata))[2]
    SegExprData <- numeric()
    for (sampleNo in 1:dim(calls(CNdata))[2]) {
        SegExpr <- numeric()
        SegData <- segmented(CNdata[, sampleNo])
        segments <- makeSegments(segmented(CNdata[, sampleNo]))
        for (j in 1:dim(segments)[1]) {
            ids <- c(segments[j, 2]:segments[j, 3])
            medSegExpr <- median(exprs(GEdata)[ids, sampleNo])
            SegExpr <- c(SegExpr, rep(medSegExpr, length(ids)))
        }
        SegExprData <- cbind(SegExprData, SegExpr)
    }
    regDetails <- makeRegions(CNprobs)
    regCalls <- numeric()
    regSegExprs <- numeric()
    regChr <- numeric()
    for (j in 1:dim(regDetails)[1]) {
        regCalls <- rbind(regCalls, calls(CNdata)[regDetails[j, 
            1], ])
        regSegExprs <- rbind(regSegExprs, SegExprData[regDetails[j, 
            1], ])
        regChr <- c(regChr, fData(CNdata)[regDetails[j, 1], 1])
    }
    chrInd <- rep(0, length(regChr))
    chrInd[(regChr%%2 == 0)] <- 1
    chrColor <- rep("blue", length(regChr))
    chrColor[(regChr%%2 == 0)] <- c("yellow")
    Y <- rep(FALSE, length(regChr))
    for (i in 2:length(regChr)) {
        if ((regChr[i - 1] != regChr[i])) {
            Y[i] <- TRUE
        }
    }
    Y[1] <- TRUE
    beginChr <- rep("", length(regChr))
    beginChr[Y] <- regChr[Y]
    CNcolor.coding <- c("red", "black", "green", "white")[1:nclass]
    def.par <- par
    fl <- layout(matrix(c(1, 2, 3, 1, 2, 3, 1, 2, 3), 3, 3, byrow = TRUE), 
        width = c(1, 9, 9))
    par(mar = c(3, 2, 4, 0))
    image(z = matrix(chrInd, nrow = 1), xaxt = "n", yaxt = "n", 
        col = c("blue", "yellow"))
    axis(2, at = (which(Y) - 1)/(length(Y) - 1), labels = regChr[Y], 
        tick = FALSE, las = 1)
    par(mar = c(3, 1, 4, 1))
    image(z = t(regCalls), xaxt = "n", yaxt = "n", col = CNcolor.coding, 
        main = "copy number data")
    par(mar = c(3, 1, 4, 1))
    image(z = t(regSegExprs), xaxt = "n", yaxt = "n", col = collist, 
        breaks = breaks, main = "gene expression data")
    par(def.par)
    return(invisible(NULL))
}

# plot heatmaps
pdf(file='matched-cn-and-expression-heatmap.pdf')
intCNGEan.heatmaps.modified(matched$CNdata.matched, matched$GEdata.matched, location='median')
dev.off()

# separate elements from the resulting object

dat3 <- data.frame(matched$CNdata.matched@featureData@data, matched$GEdata.matched@featureData@data)
colnames(dat3) <- c('chromosome', 'cn.start', 'cn.end', 'exp.probe','exp.start', 'exp.end')
dat3$exp.probe <- rownames(matched$GEdata.matched@featureData@data)

if ('symbol' %in% colnames(dat_exp))
  dat3$symbol <- dat_exp[dat3$exp.probe, 'symbol']
if ('description' %in% colnames(dat_exp))
  dat3$description <- dat_exp[dat3$exp.probe, 'description']

dat3$loss.freq <- round(mean(as.data.frame(t(assayDataElement(matched$CNdata.matched, "calls")==-1))), digits=3)
dat3$gain.freq <- round(mean(as.data.frame(t(assayDataElement(matched$CNdata.matched, "calls")==1))), digits=3)

exprs <- assayDataElement(matched$GEdata.matched, 'exprs')
samples <- microarrays <- sprintf('microarray%.3i', 1:ncol(exprs))
colnames(exprs) <- paste('chip.', samples, sep='')
dat3 <- cbind(dat3, exprs)

calls <- assayDataElement(matched$CNdata.matched, 'calls')
colnames(calls) <- paste('flag.', samples, sep='')
dat3 <- cbind(dat3, calls)

copynumber <- assayDataElement(matched$CNdata.matched, 'copynumber')
colnames(copynumber) <- paste('copynumber.', samples, sep='')
dat3 <- cbind(dat3, copynumber)

segmented <- assayDataElement(matched$CNdata.matched, 'segmented')
colnames(segmented) <- paste('segmented.', samples, sep='')
dat3 <- cbind(dat3, segmented)

probloss <- assayDataElement(matched$CNdata.matched, 'probloss')
colnames(probloss) <- paste('probloss.', samples, sep='')
dat3 <- cbind(dat3, probloss)

probnorm <- assayDataElement(matched$CNdata.matched, 'probnorm')
colnames(probnorm) <- paste('probnorm.', samples, sep='')
dat3 <- cbind(dat3, probnorm)

probgain <- assayDataElement(matched$CNdata.matched, 'probgain')
colnames(probgain) <- paste('probgain.', samples, sep='')
dat3 <- cbind(dat3, probgain)

dat3$chromosome <- as.character(dat3$chromosome)
dat3$chromosome[dat3$chromosome=='23'] <- 'X'
dat3$chromosome[dat3$chromosome=='24'] <- 'Y'
dat3$chromosome[dat3$chromosome=='25'] <- 'MT'

# add additional phenodata columns from phenodata_cgh to phenodata_exp
phenodata_cgh <- phenodata_cgh[index_cgh,]
phenodata_exp <- phenodata_exp[index_exp,]
phenodata_exp$description_cgh <- phenodata_cgh$description
for (col in setdiff(colnames(phenodata_cgh), colnames(phenodata_exp)))
  phenodata_exp[,col] <- phenodata_cgh[,col]
phenodata_exp$sample <- samples
phenodata_exp$n <- NULL

# write output
write.table(format(dat3, scientific=FALSE), file='matched-cn-and-expression.tsv', quote=FALSE, sep='\t')
write.table(format(phenodata_exp, scientific=FALSE), file='matched-phenodata.tsv', quote=FALSE, sep='\t', na='', row.names=FALSE)

# EOF
