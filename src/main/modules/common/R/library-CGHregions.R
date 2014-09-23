# Construction end deconstruction of CGHregions objects.
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-09-22

library(Biobase)
library(CGHregions)
library(WECCA)

toCghRegions <- function(df) {
  pos <- c("chromosome", "start", "end")
  if (length(setdiff(pos, colnames(df))) != 0)
    stop("CHIPSTER-NOTE: The input file is required to contain the following columns: chromosome, start, end.")
  fData <- df[, pos]
  fData$chromosome <- chromosomeToInteger(fData$chromosome)
  colnames(fData) <- c("Chromosome", "Start", "End")
  if ("num.probes" %in% colnames(df)) {
    fData$Nclone <- df$num.probes
  } else {
    fData$Nclone <- NA
  }
  fData$AveDist <- NA
  regions <- as.matrix(df[, grep("^flag\\.", colnames(df)), drop=FALSE])
  colnames(regions) <- sub("^flag\\.", "", colnames(regions))
  new("cghRegions", assayData=assayDataNew(regions=regions), featureData=new("AnnotatedDataFrame", data=fData))
}

toRegioning <- function(df) {
  pos <- c("chromosome", "start", "end")
  if (length(setdiff(pos, colnames(df))) != 0)
    stop("CHIPSTER-NOTE: The input file is required to contain the following columns: chromosome, start, end.")
  fData <- df[, pos]
  fData$chromosome <- chromosomeToInteger(fData$chromosome)
  colnames(fData) <- c("Chromosome", "Start", "End")
  if ("num.probes" %in% colnames(df)) {
    fData$Nclone <- df$num.probes
  } else {
    fData$Nclone <- NA
  }
  fData$AveDist <- NA
  hardcalls <- as.matrix(df[, grep("^flag\\.", colnames(df)), drop=FALSE])
  softcalls <- as.matrix(df[, grep("^prob", colnames(df)), drop=FALSE])
  colnames(softcalls) <- sub('\\.', '_', colnames(softcalls))
  list(ann=fData, hardcalls=hardcalls, softcalls=softcalls)
}

fromRegioning <- function(object) {
  df <- data.frame(object[["ann"]])
  colnames(df)[1:5] <- c('chromosome', 'start', 'end', 'num.probes', 'ave.dist')
  df$chromosome <- chromosomeToCharacter(df$chromosome)
  df$ave.dist <- NULL
  hardcalls <- object[["hardcalls"]]
  # calculate frequencies
  if (-2 %in% hardcalls)
    df$del.freq <- round(rowMeans(hardcalls == -2), digits=3)
  df$loss.freq <- round(rowMeans(hardcalls == -1), digits=3)
  df$gain.freq <- round(rowMeans(hardcalls == 1), digits=3)
  if (2 %in% hardcalls)
    df$amp.freq <- round(rowMeans(hardcalls == 2), digits=3)
  samples <- colnames(hardcalls)
  colnames(hardcalls) <- paste0("flag.", samples)
  df <- cbind(df, hardcalls)
  medians <- object[["medians"]]
  colnames(medians) <- paste0("chip.", samples)
  df <- cbind(df, medians)
  colnames(medians) <- paste0("segmented.", samples)
  df <- cbind(df, medians)
  softcalls <- object[["softcalls"]]
  colnames(softcalls) <- sub("_", ".", colnames(softcalls))
  df <- cbind(df, softcalls)
  df
}

# Overriding default versions.

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
  featureNames(regions) <- sprintf("%s:%i-%i", chromosomeToCharacter(chromosomes(regions)), bpstart(regions), bpend(regions))
  regions
}
environment(CGHregionsPlus) <- environment(CGHregions:::CGHregions)
CGHregions <- CGHregionsPlus

regioningPlus <- function (cghdata.called, threshold=0.00001, cghdata.regions=NULL) {
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
      } else {
        index.temp <- i + 1
        j <- j + 1
        splitter[[j]] <- index.temp
      }
    }
    region.details <- NULL
    for (i in 1:length(splitter)) {
      region.details <- rbind(region.details, c(min(splitter[[i]]), max(splitter[[i]])))
    }
    modus <- which.max(region.details[, 2] - region.details[, 1] + 1)
    return(x[region.details[modus[1], 1], ])
  }
  cat("CGHregions of hard call data...")
  if (is.null(cghdata.regions))
    cghdata.regions <- CGHregions(cghdata.called, averror = threshold)
  cat("...done", "\n")
  print(paste("threshold used:", threshold, sep = " "))
  calls.annotation <- pData(featureData(cghdata.called))
  regions.annotation <- pData(featureData(cghdata.regions))
  cat("Map regions to clones...")
  reg.to.clones <- list()
  counter <- 0
  for (chr in unique(regions.annotation[, 1])) {
    reg.ann.temp <- regions.annotation[regions.annotation[, 1] == chr, 1:4]
    for (r in 1:dim(reg.ann.temp)[1]) {
      counter <- counter + 1
      A1 <- which(calls.annotation[, 1] == chr)
      A2 <- which(calls.annotation[, 2] >= reg.ann.temp[r, 2])
      A3 <- which(calls.annotation[, 3] <= reg.ann.temp[r, 3])
      reg.to.clones[[counter]] <- intersect(intersect(A1, A2), A3)
    }
  }
  cat("...done", "\n")
  cghdata.probs <- numeric()
  for (i in 1:dim(calls(cghdata.called))[2]) {
    cghdata.probs <- cbind(cghdata.probs, cbind(probloss(cghdata.called)[, i], probnorm(cghdata.called)[, i], probgain(cghdata.called)[, i], probamp(cghdata.called)[, i]))
  }
  cat("Calculate mode soft call signature for each region...")
  cghdata.regprobs <- numeric()
  for (i in 1:length(reg.to.clones)) {
    cghdata.regprobs <- rbind(cghdata.regprobs, find.reg.modus(cghdata.probs[reg.to.clones[[i]], , drop = FALSE]))
  }
  cat("...done", "\n")
  softcalls.samplenames <- character()
  for (i in 1:dim(calls(cghdata.called))[2]) {
    if (dim(cghdata.regprobs)[2]/dim(calls(cghdata.called))[2] == 3) {
      softcalls.samplenames <- c(softcalls.samplenames, paste(c("probloss_", "probnorm_", "probgain_"), colnames(regions(cghdata.regions))[i], sep = ""))
    }
    if (dim(cghdata.regprobs)[2]/dim(calls(cghdata.called))[2] == 4) {
      softcalls.samplenames <- c(softcalls.samplenames, paste(c("probloss_", "probnorm_", "probgain_", "probamp_"), colnames(regions(cghdata.regions))[i], sep = ""))
    }
  }
  colnames(cghdata.regprobs) <- softcalls.samplenames
  rownames(cghdata.regprobs) <- rownames(regions(cghdata.regions))
  regdata <- list()
  regdata$ann <- regions.annotation
  regdata$hardcalls <- regions(cghdata.regions)
  regdata$softcalls <- cghdata.regprobs
  
  # calculate median logratios and segments
  region.medians <- regions(cghdata.regions)
  for (i in 1:nrow(region.medians)) {
    index <- chromosomes(cghdata.called) == chromosomes(cghdata.regions)[i] &
             bpstart(cghdata.called) >= bpstart(cghdata.regions)[i] &
             bpend(cghdata.called) <= bpend(cghdata.regions)[i]
    region.medians[i,] <- apply(copynumber(cghdata.called)[index, , drop=FALSE], 2, median, na.rm=TRUE)
  }
  regdata$medians <- region.medians
  
  return(regdata)
}
environment(regioningPlus) <- environment(WECCA:::regioning)
regioning <- regioningPlus

WECCA.heatmapPlus <- function (cghdata.regioned, dendrogram, build='GRCh37', ...) {
  nclass <- dim(cghdata.regioned$softcalls)[2]/dim(cghdata.regioned$hardcalls)[2]

  cols <- c('lightgreen', 'darkgreen', 'lightgray', 'darkslategray')
  chr.color <- rep(1, nrow(cghdata.regioned$ann))
  centromeres <- CGHbase:::.getCentromere(build)
  for (chr in unique(cghdata.regioned$ann$Chromosome))
    chr.color[cghdata.regioned$ann$Chromosome == chr & (cghdata.regioned$ann$Start + cghdata.regioned$ann$End) / 2 > centromeres[chr]] <- 2
  even <- cghdata.regioned$ann$Chromosome %% 2 == 0
  chr.color[even] <- chr.color[even] + 2
  chr.color <- cols[chr.color]

  Y <- rep(FALSE, dim(cghdata.regioned$hardcalls)[1])
  for (i in 2:(dim(cghdata.regioned$ann)[1])) {
    if ((cghdata.regioned$ann[i - 1, 1] != cghdata.regioned$ann[i, 1])) {
      Y[i] <- TRUE
    }
  }
  Y[1] <- TRUE
  begin.chr <- rep("", dim(cghdata.regioned$ann)[1])
  begin.chr[Y] <- cghdata.regioned$ann[Y, 1]
  color.coding <- c("red", "black", "blue", "white")[1:nclass]
  if (class(dendrogram) == "hclust")
    dendrogram <- as.dendrogram(dendrogram)
  heatmap(cghdata.regioned$hardcalls, Colv=dendrogram, Rowv=NA, col=color.coding, labRow=begin.chr, RowSideColors=chr.color, scale="none", ...)
}
environment(WECCA.heatmapPlus) <- environment(WECCA:::WECCA.heatmap)
WECCA.heatmap <- WECCA.heatmapPlus

# Additional functions.

CGHregionsManual <- function(cgh, df) {
  rownames(df) <- sprintf("%s:%i-%i", df$chromosome, df$start, df$end)
  pos <- c("chromosome", "start", "end")
  fData <- df[, pos]
  fData$chromosome <- chromosomeToInteger(fData$chromosome)
  colnames(fData) <- c("Chromosome", "Start", "End")
  fData$Nclone <- NA
  fData$AveDist <- NA
  ctdat <- data.frame(ind=1:nrow(cgh), calls(cgh))
  regcalls <- matrix(nrow=nrow(fData), ncol=ncol(cgh),
    dimnames=list(rownames(fData), sampleNames(cgh)))
  for (i in 1:nrow(fData)) {
    startend <- range(which(chromosomes(cgh) == fData[i, "Chromosome"] &
    bpstart(cgh) >= fData[i, "Start"] &
    bpend(cgh) <= fData[i, "End"]))
    fData$Nclone[i] <- startend[2] - startend[1] + 1
    signature <- CGHregions:::.whichsign2(startend, ctdat, c(-2, -1, 0, 1, 2))
    regcalls[i, ] <- as.numeric(signature[[1]][1,])
  }
  new("cghRegions", assayData=assayDataNew(regions=regcalls), featureData=new("AnnotatedDataFrame", data=fData))
}

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
  ### number of samples
  mtext(paste(ncol(x), "samples"), side=3, line=0, adj=1, cex=par("cex"))
})

# EOF
