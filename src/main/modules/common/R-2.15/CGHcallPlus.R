library(CGHcall)
library(CGHregions)
library(WECCA)

setMethod("frequencyPlot", signature(x="cghCall", y="missing"), frequencyPlotCalls)

setMethod("frequencyPlot", signature(x="cghRegions", y="missing"),
function (x, y, main='Frequency Plot', gaincol='blue', losscol='red', misscol=NA, build='GRCh37',... )
{
  chrom <- chromosomes(x)
  pos <- bpstart(x)
  pos2 <- bpend(x)
  uni.chrom <- unique(chrom)
  chrom.lengths <- CGHbase:::.getChromosomeLengths(build)[as.character(uni.chrom)]
  chrom.ends <- integer()
  cumul <- 0
  for (j in uni.chrom) {
    pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
    pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
    cumul <- cumul + chrom.lengths[as.character(j)]
    chrom.ends <- c(chrom.ends, cumul)
  }
  names(chrom.ends) <- names(chrom.lengths)
  calls <- regions(x)
  loss.freq <- rowMeans(calls < 0)
  gain.freq <- rowMeans(calls > 0)
  plot(NA, xlim=c(0, max(pos2)), ylim=c(-1,1), type='n', xlab='chromosomes', ylab='frequency', xaxs='i', xaxt='n', yaxs='i', yaxt='n', main=main,...)
  if (!is.na(misscol)) {
    rect(0, -1, max(pos2), 1, col=misscol, border=NA)
    rect(pos, -1, pos2, 1, col='white', border=NA)
  }
  rect(pos, 0, pos2, gain.freq, col=gaincol, border=gaincol)
  rect(pos, 0, pos2, -loss.freq, col=losscol, border=losscol)
  box()
  abline(h=0)
  if (length(chrom.ends) > 1)
    for (j in names(chrom.ends)[-length(chrom.ends)])
      abline(v=chrom.ends[j], lty='dashed')
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1,at=ax,labels=uni.chrom,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
  mtext('gains', side=2, line=3, at=0.5)
  mtext('losses', side=2, line=3, at=-0.5)
  ### number of data points
  mtext(paste(nrow(x), 'regions'), side=3, line=0, adj=0)
})

segmentDataWeighted <-
function (input, method = "DNAcopy", ...) 
{
    if (method == "DNAcopy") {
        CNA.object <- DNAcopy::CNA(copynumber(input), chromosomes(input), 
            bpstart(input), data.type = "logratio")
        cat("Start data segmentation .. \n")
        segmented <- segment(CNA.object, ...)
        numclone <- segmented$output$num.mark
        smrat <- segmented$output$seg
        numsmrat <- cbind(smrat, numclone)
        repdata <- function(row) {
            rep(row[1], row[2])
        }
        makelist <- apply(numsmrat, 1, repdata)
        joined <- unlist(makelist)
        rm(makelist)
        joined <- matrix(joined, ncol = ncol(input), byrow = FALSE)
        joined <- CGHcall:::.assignNames(joined, input)
        result <- CGHcall:::.segFromRaw(input, joined)
    }
    result
}

CGHregionsPlus <- function(input, ...) {
  regions <- CGHregions:::CGHregions(input, ...)
  # End positions of regions should be the end position of the last data point of that region,
  # but instead CGHregions returns the start position of the last data point.
  # Check if that is indeed the case:
  if (class(input) == 'cghCall') {
    if (sum(regions@featureData@data$End %in% input@featureData@data$Start) > sum(regions@featureData@data$End %in% input@featureData@data$End))
      for (row in rownames(regions@featureData@data))
        regions@featureData@data[row, 'End'] <- input@featureData@data[input@featureData@data$Chromosome == regions@featureData@data[row, 'Chromosome'] & input@featureData@data$Start == regions@featureData@data[row, 'End'], 'End'][1]
  }
  regions
}
environment(CGHregionsPlus) <- environment(CGHregions:::CGHregions)
CGHregions <- CGHregionsPlus

regioningPlus <- function (cghdata.called, threshold = 0.00001, cghdata.regions = NULL) 
{
    find.reg.modus <- function(x) {
        if (nrow(x) == 1)
            return(x)
        splitter <- list()
        splitter[[1]] <- c(1)
        index.temp <- 1
        j <- 1
        for (i in 1:(dim(x)[1] - 1)) {
            if (all(x[i, ] == x[i + 1, ])) {
                index.temp <- c(index.temp, i + 1)
                splitter[[j]] <- index.temp
            }
            else {
                index.temp <- i + 1
                j <- j + 1
                splitter[[j]] <- index.temp
            }
        }
        region.details <- NULL
        for (i in 1:length(splitter)) {
            region.details <- rbind(region.details, c(min(splitter[[i]]), 
                max(splitter[[i]])))
        }
        modus <- which.max(region.details[, 2] - region.details[, 
            1] + 1)
        return(x[region.details[modus[1], 1], ])
    }
    cat("CGHregions of hard call data...")
    if (is.null(cghdata.regions))
        cghdata.regions <- CGHregionsPlus(cghdata.called, averror = threshold)
    cat("...done", "\n")
    print(paste("threshold used:", threshold, sep = " "))
    calls.annotation <- pData(featureData(cghdata.called))
    regions.annotation <- pData(featureData(cghdata.regions))
    cat("Map regions to clones...")
    reg.to.clones <- list()
    counter <- 0
    for (chr in unique(regions.annotation[, 1])) {
        reg.ann.temp <- regions.annotation[regions.annotation[, 
            1] == chr, 1:4]
        for (r in 1:dim(reg.ann.temp)[1]) {
            counter <- counter + 1
            A1 <- which(calls.annotation[, 1] == chr)
            A2 <- which(calls.annotation[, 2] >= reg.ann.temp[r, 
                2])
            A3 <- which(calls.annotation[, 3] <= reg.ann.temp[r, 
                3])
            reg.to.clones[[counter]] <- intersect(intersect(A1, 
                A2), A3)
        }
    }
    cat("...done", "\n")
    cghdata.probs <- numeric()
    for (i in 1:dim(calls(cghdata.called))[2]) {
        cghdata.probs <- cbind(cghdata.probs, cbind(probloss(cghdata.called)[, 
            i], probnorm(cghdata.called)[, i], probgain(cghdata.called)[, 
            i], probamp(cghdata.called)[, i]))
    }
    cat("Calculate mode soft call signature for each region...")
    cghdata.regprobs <- numeric()
    for (i in 1:length(reg.to.clones)) {
        cghdata.regprobs <- rbind(cghdata.regprobs, find.reg.modus(cghdata.probs[reg.to.clones[[i]], 
            , drop = FALSE]))
    }
    cat("...done", "\n")
    softcalls.samplenames <- character()
    for (i in 1:dim(calls(cghdata.called))[2]) {
        if (dim(cghdata.regprobs)[2]/dim(calls(cghdata.called))[2] == 
            3) {
            softcalls.samplenames <- c(softcalls.samplenames, 
                paste(c("probloss_", "probnorm_", "probgain_"), 
                  colnames(regions(cghdata.regions))[i], sep = ""))
        }
        if (dim(cghdata.regprobs)[2]/dim(calls(cghdata.called))[2] == 
            4) {
            softcalls.samplenames <- c(softcalls.samplenames, 
                paste(c("probloss_", "probnorm_", "probgain_", 
                  "probamp_"), colnames(regions(cghdata.regions))[i], 
                  sep = ""))
        }
    }
    colnames(cghdata.regprobs) <- softcalls.samplenames
    rownames(cghdata.regprobs) <- rownames(regions(cghdata.regions))
    regdata <- list()
    regdata$ann <- regions.annotation
    regdata$hardcalls <- regions(cghdata.regions)
    regdata$softcalls <- cghdata.regprobs
    return(regdata)
}
environment(regioningPlus) <- environment(WECCA:::regioning)
regioning <- regioningPlus

WECCA.heatmapPlus <- function (cghdata.regioned, dendrogram, build='GRCh37', ...) 
{
    nclass <- dim(cghdata.regioned$softcalls)[2]/dim(cghdata.regioned$hardcalls)[2]

    cols <- c('lightgreen', 'darkgreen', 'lightgray', 'darkslategray')
    chr.color <- rep(1, nrow(cghdata.regioned$ann))
    centromeres <- CGHbase:::.getCentromere(build)
    for (chr in unique(cghdata.regioned$ann$Chromosome))
      chr.color[cghdata.regioned$ann$Chromosome == chr &
               (cghdata.regioned$ann$Start + cghdata.regioned$ann$End) / 2 > centromeres[chr]] <- 2
    even <- cghdata.regioned$ann$Chromosome %% 2 == 0
    chr.color[even] <- chr.color[even] + 2
    chr.color <- cols[chr.color]

    Y <- rep(FALSE, dim(cghdata.regioned$hardcalls)[1])
    for (i in 2:(dim(cghdata.regioned$ann)[1])) {
        if ((cghdata.regioned$ann[i - 1, 1] != cghdata.regioned$ann[i, 
            1])) {
            Y[i] <- TRUE
        }
    }
    Y[1] <- TRUE
    begin.chr <- rep("", dim(cghdata.regioned$ann)[1])
    begin.chr[Y] <- cghdata.regioned$ann[Y, 1]
    color.coding <- c("red", "black", "blue", "white")[1:nclass]
    heatmap(cghdata.regioned$hardcalls, Colv = as.dendrogram(dendrogram), 
        Rowv = NA, col = color.coding, labRow = begin.chr, RowSideColors = chr.color, 
        scale = "none", ...)
}
environment(WECCA.heatmapPlus) <- environment(WECCA:::WECCA.heatmap)
WECCA.heatmap <- WECCA.heatmapPlus

# EOF
