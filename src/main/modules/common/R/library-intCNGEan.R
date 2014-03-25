# Construction end deconstruction of ExpressionSet objects.
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-25

library(Biobase)
library(intCNGEan)

toIntCNGEanMatched <- function(df) {
  # check that the input file seems to be coming from the script used to match the two data sets
  pos <- c("chromosome","cn.start","cn.end","exp.start","exp.end")
  if (length(setdiff(pos, colnames(df))) != 0)
    stop("CHIPSTER-NOTE: This tool can only be run on the output file from the tool Match copy number and expression probes (matched-cn-and-expression.tsv).")

  df$chromosome <- chromosomeToInteger(df$chromosome)
  
  exprs <- as.matrix(df[, grep("^chip\\.", colnames(df))])
  samples <- sub("^chip\\.", "", colnames(exprs))
  colnames(exprs) <- samples
  rownames(exprs) <- df$exp.probe
  calls <- as.matrix(df[, grep("^flag\\.", colnames(df))])
  colnames(calls) <- samples
  copynumber <- as.matrix(df[, grep("^logratio\\.", colnames(df))])
  colnames(copynumber) <- samples
  segmented <- as.matrix(df[, grep("^segmented\\.", colnames(df))])
  colnames(segmented) <- samples
  probloss <- as.matrix(df[, grep("^probloss\\.", colnames(df))])
  colnames(probloss) <- samples
  probnorm <- as.matrix(df[, grep("^probnorm\\.", colnames(df))])
  colnames(probnorm) <- samples
  probgain <- as.matrix(df[, grep("^probgain\\.", colnames(df))])
  colnames(probgain) <- samples

  cgh <- new("cghCall", assayData=assayDataNew(calls=calls, copynumber=copynumber, segmented=segmented, probloss=probloss, probnorm=probnorm, probgain=probgain), featureData=new("AnnotatedDataFrame", data=data.frame(Chromosome=df$chromosome, Start=df$cn.start, End=df$cn.end, row.names=row.names(df))))
  exp <- new("ExpressionSet", exprs=exprs, featureData=new("AnnotatedDataFrame", data=data.frame(Chromosome=df$chromosome, Start=df$exp.start, End=df$exp.end, row.names=df$exp.probe)))

  list(CNdata.matched=cgh, GEdata.matched=exp)
}

toIntCNGEanTuned <- function(df) {
  df <- df[order(df$gene.id),]
  datafortest <- as.matrix(df[, grep("^(prob1|prob2|chip)\\.", colnames(df))])
  nosamp <- ncol(datafortest)/3
  rownames(datafortest) <- df$probes
  colnames(datafortest) <- sub("^chip\\.", "", colnames(datafortest))
  colnames(datafortest)[1:(2*nosamp)] <- ""
  callprobs <- as.matrix(df[, c("av.probs.1", "av.probs.2")])
  rownames(callprobs) <- df$probes
  colnames(callprobs) <- NULL
  lossorgain <- df$comparison
  names(lossorgain) <- df$probes
  genestotest <- df$gene.id
  names(genestotest) <- df$probes
  ann <- df[, c("chromosome", "start", "end")]
  colnames(ann) <- c("Chromosome", "Start", "End")
  alleffects <- df$effect.size
  tuned <- list(datafortest=datafortest, lossorgain=lossorgain, genestotest=genestotest, callprobs=callprobs, alleffects=alleffects, ann=ann, nosamp=nosamp)
}

fromIntCNGEanMatched <- function(object) {
  df1 <- data.frame(fData(object$CNdata.matched))
  colnames(df1)[1:3] <- c("chromosome", "cn.start", "cn.end")
  df2 <- data.frame(fData(object$GEdata.matched))
  colnames(df2)[1:3] <- c("exp.probe", "exp.start", "exp.end")
  df2$exp.probe <- featureNames(object$GEdata.matched)
  df <- cbind(df1, df2)

  calls <- assayDataElement(object$CNdata.matched, "calls")
  df$loss.freq <- round(rowMeans(calls == -1), digits=3)
  df$gain.freq <- round(rowMeans(calls == 1), digits=3)

  exprs <- assayDataElement(object$GEdata.matched, "exprs")
  samples <- microarrays <- sprintf("microarray%.3i", 1:ncol(exprs))
  colnames(exprs) <- paste("chip.", samples, sep="")

  colnames(calls) <- paste("flag.", samples, sep="")

  copynumber <- assayDataElement(object$CNdata.matched, "copynumber")
  colnames(copynumber) <- paste("logratio.", samples, sep="")

  segmented <- assayDataElement(object$CNdata.matched, "segmented")
  colnames(segmented) <- paste("segmented.", samples, sep="")

  probloss <- assayDataElement(object$CNdata.matched, "probloss")
  colnames(probloss) <- paste("probloss.", samples, sep="")

  probnorm <- assayDataElement(object$CNdata.matched, "probnorm")
  colnames(probnorm) <- paste("probnorm.", samples, sep="")

  probgain <- assayDataElement(object$CNdata.matched, "probgain")
  colnames(probgain) <- paste("probgain.", samples, sep="")

  df <- cbind(df, exprs, calls, copynumber, segmented, probloss, probnorm, probgain)

  df$chromosome <- chromosomeToCharacter(df$chromosome)
  df
}

intCNGEan.heatmapsPlus <- function (CNdata, GEdata, location = "mode", colorbreaks = "equiquantiles") 
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
environment(intCNGEan.heatmapsPlus) <- environment(intCNGEan.heatmaps)
intCNGEan.heatmaps <- intCNGEan.heatmapsPlus

intCNGEan.tunePlus <- function (CNdata, GEdata, test.statistic, ngenetune = 250, nperm_tuning = 250, 
    minCallProbMass = 0.05, trace = TRUE) 
{
    countth <- function(statlist) {
        threshold <- as.numeric(statlist[length(statlist)])
        statlist <- as.numeric(statlist[c(1:(length(statlist) - 
            1))])
        return(length(statlist[statlist >= threshold]))
    }
    pval.perm.marg <- function(observed, permuted, nperm) {
        perm.and.obs <- cbind(permuted, observed)
        return(apply(perm.and.obs, 1, "countth")/nperm)
    }
    rawps <- function(stats.obs, nulldists, nperm) {
        pval.ln <- pval.perm.marg(stats.obs, nulldists, nperm)
        return(pval.ln)
    }
    wcvm.test.stats <- function(cgh.em, nosamp, a) {
        cgh.2cat <- matrix(cgh.em[c(1:(a * nosamp))], ncol = a, 
            byrow = TRUE)
        alphaas <- t(cgh.2cat) %*% cgh.2cat
        cs <- as.numeric(solve(alphaas) %*% matrix(c(-1, 1), 
            ncol = 1))
        cgh.em <- cbind(cgh.2cat, cgh.em[c((a * nosamp + 1):((a + 
            1) * nosamp))])
        cgh.em <- cbind(cgh.em[order(cgh.em[, 3]), ], rep(1/dim(cgh.em)[1], 
            dim(cgh.em)[1]))
        cgh.em <- cbind(cgh.em, cumsum(cgh.em[, 4]), cumsum(cgh.em[, 
            1] * cgh.em[, dim(cgh.em)[2]]) * cs[1], cumsum(cgh.em[, 
            2] * cgh.em[, dim(cgh.em)[2]]) * cs[2])
        test.stat <- -sum(cgh.em[, 7] + cgh.em[, 6])/nosamp
        return(test.stat)
    }
    wmw.test.stats <- function(cgh.em, nosamp, a) {
        cgh.em <- cbind(matrix(cgh.em[c(1:(a * nosamp))], ncol = a, 
            byrow = TRUE), cgh.em[c((a * nosamp + 1):((a + 1) * 
            nosamp))])
        cgh.em <- cgh.em[order(cgh.em[, (a + 1)]), ]
        cgh.em <- cbind(cgh.em, rbind(rep(0, a), apply(cgh.em[, 
            1:a], 2, cumsum)[-nosamp, ]))
        test.stat <- sum(cgh.em[, a + 2] * cgh.em[, 2])
        return(test.stat)
    }
    prob.test.stats <- function(cgh.em, nosamp, a) {
        cgh.2cat <- matrix(cgh.em[c(1:(a * nosamp))], ncol = a, 
            byrow = TRUE)
        alphaas <- t(cgh.2cat) %*% cgh.2cat
        cs <- c(det(alphaas), alphaas[1, 2] * sum(alphaas)/2)
        cgh.em <- cbind(matrix(cgh.em[c(1:(a * nosamp))], ncol = a, 
            byrow = TRUE), cgh.em[c((a * nosamp + 1):((a + 1) * 
            nosamp))])
        cgh.em <- cgh.em[order(cgh.em[, (a + 1)]), ]
        cgh.em <- cbind(cgh.em, rbind(rep(0, a), apply(cgh.em[, 
            1:a], 2, cumsum)[-nosamp, ]))
        test.stat <- (sum(cgh.em[, a + 2] * cgh.em[, 2]) - cs[2])/cs[1]
        return(test.stat)
    }
    pretest <- function(alphascgh) {
        alphas <- alphascgh[1:3]
        probs <- matrix(alphascgh[-(1:3)], ncol = 3, byrow = TRUE)
        if (alphas[1] >= alphas[3]) {
            probs2 <- cbind(probs[, 1], probs[, 2] + probs[, 
                3])
            alphas2 <- c(alphas[1], alphas[2] + alphas[3])
            return(c(1, alphas2, as.vector(t(probs2))))
        }
        else {
            probs2 <- cbind(probs[, 1] + probs[, 2], probs[, 
                3])
            alphas2 <- c(alphas[1] + alphas[2], alphas[3])
            return(c(2, alphas2, as.vector(t(probs2))))
        }
    }
    alphaest <- function(cgh.em, nosamp, a) {
        cgh.em <- matrix(cgh.em[c(1:(a * nosamp))], ncol = a, 
            byrow = TRUE)
        return(apply(cgh.em, 2, mean))
    }
    alphabivariate <- function(cgh.em, nosamp, a) {
        cgh.em <- matrix(cgh.em[c(1:(a * nosamp))], ncol = a, 
            byrow = TRUE)
        cgh.em1 <- cgh.em[, 1]
        cgh.em2 <- cgh.em[, 2]
        return(c(1/nosamp * (cgh.em1 %*% cgh.em1), 1/nosamp * 
            (cgh.em2 %*% cgh.em2), 1/nosamp * (cgh.em1 %*% cgh.em2)))
    }
    probs2calls <- function(problist) {
        maxprobposition <- which.max(problist)
        call.list <- rep(0, length(problist))
        call.list[maxprobposition] <- 1
        return(call.list)
    }
    shift.est <- function(row, nosamp, a, data2, alphabmat, alphanmat, 
        minalphathr) {
        cgh.em <- cbind(matrix(data2[row, c(1:(a * nosamp))], 
            ncol = a, byrow = TRUE), data2[row, c((a * nosamp + 
            1):((a + 1) * nosamp))])
        alphasbiv <- alphabmat[row, ]
        alphas <- alphanmat[row, ]
        c1 <- (alphasbiv[1]/alphasbiv[3] - alphas[1]/alphas[2])^(-1)
        c2 <- (alphasbiv[2]/alphasbiv[3] - alphas[2]/alphas[1])^(-1)
        if (is.na(c1) | is.na(c2)) {
            return(NA)
        }
        else {
            if (max(c1, c2) >= minalphathr) {
                return(NA)
            }
            else {
                mu1 <- 1/nosamp * sum(cgh.em[, 3] * (cgh.em[, 
                  1] * alphas[2] - alphasbiv[3])/(alphasbiv[1] * 
                  alphas[2] - alphas[1] * alphasbiv[3]))
                mu2 <- 1/nosamp * sum(cgh.em[, 3] * (cgh.em[, 
                  2] * alphas[1] - alphasbiv[3])/(alphasbiv[2] * 
                  alphas[1] - alphas[2] * alphasbiv[3]))
                shiftest <- mu2 - mu1
                return(shiftest)
            }
        }
    }
    powerunbalance <- function(row, alphabmat) {
        alphasbiv <- alphabmat[row, ]
        unbalance <- alphasbiv[1] * alphasbiv[2] - alphasbiv[3]^2
        return(unbalance)
    }
    countdiscoveries <- function(powerquant, rawp_powerunbal, 
        fdrcutoff) {
        rawp_filt <- rawp_powerunbal[rawp_powerunbal[, 2] >= 
            powerquant, ][, 1]
        m <- length(rawp_filt)
        adpjrawp_filt <- cbind(rawp_filt, p.adjust(rawp_filt, 
            "BH"))
        selected <- matrix(adpjrawp_filt[adpjrawp_filt[, 2] <= 
            fdrcutoff, ], ncol = 2)
        count <- nrow(selected)
        false <- ifelse(count > 0, max(selected[, 1]) * m, 0)
        countS <- count - false
        return(countS)
    }
    redraw <- function(col, shiftest, shiftsam, alpharow, dataexp, 
        datacgh) {
        if (is.na(shiftest)) {
            return(sample(dataexp, 1))
        }
        else {
            alpha <- alpharow[1]
            Ig <- sample(c(1, 2), 1, prob = c(alpha, 1 - alpha))
            if (Ig == 1) {
                xib <- sample(dataexp, 1, prob = datacgh[, 1]) + 
                  shiftest
            }
            else {
                xib <- sample(dataexp, 1, prob = datacgh[, 2])
            }
            Ig2 <- sample(c(1, 2), 1, prob = c(max(0, min(1, 
                datacgh[col, 1])), max(0, min(1, 1 - datacgh[col, 
                1]))))
            if (Ig2 == 1) {
                xib2 <- xib
            }
            else {
                xib2 <- xib + shiftsam
            }
            return(xib2)
        }
    }
    redrawcol <- function(row, allest, alphanmat, shiftsampled, 
        data2, a, nosamp) {
        shiftest <- allest[row]
        datacgh <- matrix(data2[row, c(1:(a * nosamp))], ncol = a, 
            byrow = TRUE)
        dataexp <- data2[row, c((a * nosamp + 1):((a + 1) * nosamp))]
        shiftsam <- shiftsampled[row]
        alpharow <- alphanmat[row, ]
        result <- sapply(1:nosamp, redraw, shiftest = shiftest, 
            alpharow = alpharow, shiftsam = shiftsam, dataexp = dataexp, 
            datacgh = datacgh)
        return(result)
    }
    pi0est <- function(data.in, test.stat = "wmw", seqg = seqgenes, 
        nperm = nperm_tuning, a, nosamp) {
        data.in2 <- data.in[seqg, ]
        if (test.stat == "wmw") 
            test.stats <- apply(data.in2, 1, wmw.test.stats, 
                nosamp, a)
        else test.stats <- apply(data.in2, 1, wcvm.test.stats, 
            nosamp, a)
        rawpvals.test <- rawps(test.stats, null.dists.tune, nperm)
        pi0 <- convest(rawpvals.test)
        return(pi0)
    }
    tuning <- function(data.in, datacgh.in, test.stat = "wmw", 
        seqg = seqgenes, nperm = nperm_tuning, pi0, fdrcut = 0.25, 
        nresamp = 100, gridnr = 30, minim = 10, a = a, nosamp = nosamp, 
        shiftsam = shiftsam, trace) {
        probs <- seq(0, gridnr/(gridnr + 1), 1/(gridnr + 1))
        probstruncate <- probs[probs <= ((100 - minim)/100)]
        allest2 <- allest[seqg]
        alphanewmat2 <- alphanewmat[seqg, ]
        data.in2 = data.in[seqg, ]
        datacgh.in2 <- datacgh.in[seqg, ]
        powerunbal2 <- powerunbal[seqg]
        quant <- quantile(powerunbal2, probs = probstruncate)
        whatthrmat <- c()
        for (j in 1:nresamp) {
            if ((j%%50) == 0) {
                if (trace) {
                  cat(paste(j, "of", nresamp, "resamples done, and counting...", 
                    sep = " "), "\n")
                }
            }
            shiftsampled <- sample(shiftsam, nrow(data.in2), 
                replace = TRUE)
            pi0cutoff <- quantile(shiftsampled, pi0)
            shiftsampled <- sapply(shiftsampled, function(x) {
                if (x <= pi0cutoff) 
                  return(0)
                else return(x)
            })
            redrawall <- t(sapply(1:nrow(data.in2), redrawcol, 
                allest = allest2, alphanmat = alphanewmat2, shiftsampled = shiftsampled, 
                data2 = data.in2, a = a, nosamp = nosamp))
            newdata <- cbind(datacgh.in2[, -(1:3)], redrawall)
            if (test.stat == "wmw") 
                test.resamp <- apply(newdata, 1, wmw.test.stats, 
                  nosamp, a)
            else test.resamp <- apply(newdata, 1, wcvm.test.stats, 
                nosamp, a)
            rawpvals.test.resamp <- rawps(test.resamp, null.dists.tune, 
                nperm)
            rawp_powerunbal <- cbind(rawpvals.test.resamp, powerunbal2)
            whatthresh <- sapply(quant, countdiscoveries, rawp_powerunbal = rawp_powerunbal, 
                fdrcutoff = fdrcut)
            whatthrmat <- rbind(whatthrmat, whatthresh)
        }
        whatthrmean <- apply(whatthrmat, 2, mean)
        bestthr <- which.max(whatthrmean)
        return(quant[bestthr])
    }
    datareduce <- function(powerunbal, unbalthr) {
        datarows <- which(powerunbal >= unbalthr)
        return(datarows)
    }
    nulldist.all.wcvm <- function(data.both, nosamp, a, nperm, 
        trace) {
        cvm.like.mat <- NULL
        for (i in 1:nperm) {
            if ((i%%50) == 0) {
                if (trace) {
                  cat(i, " of ", nperm, " permutations done, and counting...", 
                    "\n")
                }
            }
            x <- sample(1:nosamp, nosamp) + a * nosamp
            data.ran <- cbind(data.both[, c(1:(a * nosamp))], 
                data.both[, x])
            cvm.like.ran <- apply(data.ran, 1, wcvm.test.stats, 
                nosamp, a)
            cvm.like.mat <- cbind(cvm.like.mat, cvm.like.ran)
        }
        return(cvm.like.mat)
    }
    nulldist.all.wmw <- function(data.both, nosamp, a, nperm, 
        trace) {
        wmw.ln.mat <- NULL
        cghdata <- data.both[, c(1:(a * nosamp))]
        for (i in 1:nperm) {
            if ((i%%50) == 0) {
                if (trace) {
                  cat(i, " of ", nperm, " permutations done, and counting...", 
                    "\n")
                }
            }
            x <- sample(1:nosamp, nosamp) + a * nosamp
            data.ran <- cbind(cghdata, data.both[, x])
            wmw.ran <- apply(data.ran, 1, wmw.test.stats, nosamp, 
                a)
            wmw.ln.mat <- cbind(wmw.ln.mat, wmw.ran)
        }
        return(wmw.ln.mat)
    }
    if (sum(is.na(exprs(GEdata))) > 0) {
        stop("Gene expression matrix contains missing values: not allowed.")
    }
    if (dim(exprs(GEdata))[2] != dim(copynumber(CNdata))[2]) {
        stop("Gene expression and copy number matrices contain unequal number of samples.")
    }
    if (dim(exprs(GEdata))[1] != dim(copynumber(CNdata))[1]) {
        stop("Gene expression and copy number matrices contain unequal number of samples: impossible after matching.")
    }
    data <- list()
    data$ann <- fData(GEdata)
    data$em <- exprs(GEdata)
    nosamp <- ncol(exprs(GEdata))
    cghdata.probs <- numeric()
    for (i in 1:dim(calls(CNdata))[2]) {
        if (is.null(probamp(CNdata)[, i])) {
            cghdata.probs <- cbind(cghdata.probs, cbind(probloss(CNdata)[, 
                i], probnorm(CNdata)[, i], probgain(CNdata)[, 
                i]))
        }
        else {
            cghdata.probs <- cbind(cghdata.probs, cbind(probloss(CNdata)[, 
                i], probnorm(CNdata)[, i], probgain(CNdata)[, 
                i] + probamp(CNdata)[, i]))
        }
    }
    data$cgh <- cghdata.probs
    nclass <- dim(cghdata.probs)[2]/dim(calls(CNdata))[2]
    if (trace) {
        cat("pre-test...", "\n")
    }
    alphamat <- t(apply(data$cgh, 1, alphaest, nosamp = nosamp, 
        a = nclass))
    datacgh2 <- as.matrix(t(apply(cbind(alphamat, data$cgh), 
        1, pretest)))
    lossorgain <- datacgh2[, 1]
    alphas2 <- datacgh2[, 2:3]
    data.both <- as.matrix(cbind(datacgh2, data$em)[, -(1:3)])
    suffCPM <- which(apply(alphas2, 1, min) > minCallProbMass)
    data$ann <- data$ann[suffCPM, ]
    data$cgh <- data$cgh[suffCPM, ]
    data$em <- data$em[suffCPM, ]
    alphamat <- alphamat[suffCPM, ]
    datacgh2 <- datacgh2[suffCPM, ]
    lossorgain <- lossorgain[suffCPM]
    alphas2 <- alphas2[suffCPM, ]
    data.both <- data.both[suffCPM, ]
    init.prop.kept <- length(suffCPM)/dim(exprs(GEdata))[1]
    alphanewmat <- alphas2
    alphabivmat <- t(apply(data.both, 1, alphabivariate, nosamp = nosamp, 
        a = 2))
    colnames(alphabivmat) <- c("a11", "a22", "a12")
    if (trace) {
        cat("tuning started...", "\n")
    }
    powerunbal <- sapply(1:nrow(data.both), powerunbalance, alphabmat = alphabivmat)
    allest <- sapply(1:nrow(data.both), shift.est, data2 = data.both, 
        alphabmat = alphabivmat, alphanmat = alphanewmat, nosamp = nosamp, 
        a = 2, minalphathr = 10)
    estrobust <- sapply(1:nrow(data.both), shift.est, data2 = data.both, 
        alphabmat = alphabivmat, alphanmat = alphanewmat, nosamp = nosamp, 
        a = 2, minalphathr = 3)
    shiftsam <- estrobust[!is.na(estrobust)]
    seqgenes <- floor(seq(1, nrow(data.both), length.out = ngenetune))
    if (test.statistic == "wcvm") {
        null.dists.tune <- nulldist.all.wcvm(data.both[seqgenes, 
            ], nosamp, 2, nperm = nperm_tuning, trace)
    }
    else {
        null.dists.tune <- nulldist.all.wmw(data.both[seqgenes, 
            ], nosamp, 2, nperm = nperm_tuning, trace)
    }
    pi0 <- pi0est(data.in = data.both, seqg = seqgenes, test.stat = test.statistic, 
        nperm = nperm_tuning, a = 2, nosamp = nosamp)
    data.tuned <- tuning(data.in = data.both, datacgh.in = datacgh2, 
        seqg = seqgenes, test.stat = test.statistic, nperm = nperm_tuning, 
        pi0 = pi0, fdrcut = 0.05, nresamp = 100, gridnr = 200, 
        minim = 10, a = 2, nosamp = nosamp, shiftsam = shiftsam, 
        trace)
    genestotest <- datareduce(powerunbal, as.double(data.tuned)) # as.real -> as.double
    prop_kept <- length(genestotest)/nrow(data.both)
    if (trace) {
        cat("ready: tuning done", "\n")
    }
    datacgh2 <- datacgh2[genestotest, ]
    lossorgain <- lossorgain[genestotest]
    alphas2 <- alphas2[genestotest, ]
    data.both <- data.both[genestotest, ]
    alphanewmat <- alphanewmat[genestotest, ]
    alphabivmat <- alphabivmat[genestotest, ]
    powerunbal <- powerunbal[genestotest]
    allest <- allest[genestotest]
    estrobust <- estrobust[genestotest]
    data.tuned <- list()
    data.tuned$datafortest <- data.both
    data.tuned$lossorgain <- lossorgain
    data.tuned$genestotest <- suffCPM[genestotest]
    data.tuned$callprobs <- alphas2
    data.tuned$alleffects <- allest
    data.tuned$ann <- data$ann[genestotest, ]
    data.tuned$nosamp <- nosamp
    return(data.tuned)
}
environment(intCNGEan.tunePlus) <- environment(intCNGEan.tune)
intCNGEan.tune <- intCNGEan.tunePlus

# EOF
