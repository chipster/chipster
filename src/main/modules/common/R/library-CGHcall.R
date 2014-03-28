# Construction end deconstruction of CGHcall objects.
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-25

library(Biobase)
library(CGHcall)
library(QDNAseq)

toCgh <- function(df, level="any", maxStates=5) {
  pos <- c("chromosome", "start", "end")
  if (length(setdiff(pos, colnames(df))) != 0)
    stop("CHIPSTER-NOTE: The input file is required to contain the following columns: chromosome, start, end.")
  df$chromosome <- chromosomeToInteger(df$chromosome)
  df <- df[!is.na(df$chromosome) & !is.na(df$start), ]
  df <- df[order(df$chromosome, df$start), ]
  fData <- df[, pos]
  colnames(fData) <- c("Chromosome", "Start", "End")

  copynumber <- as.matrix(df[, grep("^chip\\.", colnames(df)), drop=FALSE])
  samples <- sub("^chip\\.", "", colnames(copynumber))
  colnames(copynumber) <- samples
  if (level == "copynumber")
    return(new("cghRaw", copynumber=copynumber, featureData=new("AnnotatedDataFrame", data=fData)))
  segmented <- as.matrix(df[, grep("^segmented\\.", colnames(df)), drop=FALSE])
  if (ncol(segmented) == 0) {
    if (level == "any") {
      return(new("cghRaw", copynumber=copynumber, featureData=new("AnnotatedDataFrame", data=fData)))
    } else {
      stop("CHIPSTER-NOTE: The input file does not contain segmented data. Please segment it first.")
    }
  }
  colnames(segmented) <- samples
  if (level == "segmented")
    return(new("cghSeg", copynumber=copynumber, segmented=segmented, featureData=new("AnnotatedDataFrame", data=fData)))
  calls <- as.matrix(df[, grep("^flag\\.", colnames(df)), drop=FALSE])
  if (ncol(calls) == 0) {
    if (level == "any") {
      return(new("cghSeg", copynumber=copynumber, segmented=segmented, featureData=new("AnnotatedDataFrame", data=fData)))
    } else {
      stop("CHIPSTER-NOTE: The input file does not contain called data. Please call aberrations first.")
    }
  }
  colnames(calls) <- samples
  probloss <- as.matrix(df[, grep("^probloss\\.", colnames(df)), drop=FALSE])
  if (ncol(probloss) == 0)
    probloss <- (calls == -1) + 0
  colnames(probloss) <- samples
  probnorm <- as.matrix(df[, grep("^probnorm\\.", colnames(df)), drop=FALSE])
  if (ncol(probnorm) == 0)
    probnorm <- (calls == 0) + 0
  colnames(probnorm) <- samples
  probgain <- as.matrix(df[, grep("^probgain\\.", colnames(df)), drop=FALSE])
  if (ncol(probgain) == 0)
    probgain <- (calls == 1) + 0
  colnames(probgain) <- samples
  probdel <- as.matrix(df[, grep("^probdel\\.", colnames(df)), drop=FALSE])
  if (ncol(probdel) == 0 && -2 %in% calls)
    probdel <- (calls == -2) + 0
  if (maxStates < 5) {
    calls[calls == -2] <- -1
    if (ncol(probdel) > 0)
      probloss <- probloss + probdel
  }
  probamp <- as.matrix(df[, grep("^probamp\\.", colnames(df)), drop=FALSE])
  if (ncol(probamp) == 0 && 2 %in% calls)
    probamp <- (calls == 2) + 0
  if (maxStates < 4) {
    calls[calls == 2] <- 1
    if (ncol(probamp) > 0)
      probgain <- probgain + probamp
  }
  cgh <- new("cghCall", copynumber=copynumber, segmented=segmented, calls=calls, probloss=probloss, probnorm=probnorm, probgain=probgain, featureData=new("AnnotatedDataFrame", data=fData))
  if (ncol(probdel) > 0 && maxStates >= 5) {
    colnames(probdel) <- samples
    assayDataElement(cgh, "probdloss") <- probdel
  }
  if (ncol(probamp) > 0 && maxStates >= 4) {
    colnames(probamp) <- samples
    assayDataElement(cgh, "probamp") <- probamp
  }
  cgh
}

fromCgh <- function(object) {
  samples <- sampleNames(object)
  df <- fData(object)[, c("Chromosome", "Start", "End")]
  colnames(df) <- c("chromosome", "start", "end")
  if ("calls" %in% assayDataElementNames(object)) {
    calls <- assayDataElement(object, "calls")
    colnames(calls) <- paste0("flag.", colnames(calls))
    if (-2 %in% calls)
      df$del.freq <- round(rowMeans(calls == -2), digits=3)
    df$loss.freq <- round(rowMeans(calls == -1), digits=3)
    df$gain.freq <- round(rowMeans(calls == 1), digits=3)
    if (2 %in% calls)
      df$amp.freq <- round(rowMeans(calls == 2), digits=3)
    df <- cbind(df, calls)
  }
  copynumber <- signif(assayDataElement(object, "copynumber"), digits=3)
  colnames(copynumber) <- paste0("chip.", colnames(copynumber))
  df <- cbind(df, copynumber)
  if ("segmented" %in% assayDataElementNames(object)) {
    segmented <- signif(assayDataElement(object, "segmented"), digits=3)
    colnames(segmented) <- paste0("segmented.", colnames(segmented))
    df <- cbind(df, segmented)
  }
  if ("calls" %in% assayDataElementNames(object)) {
    if ("probdloss" %in% assayDataElementNames(object)) {
      probdel <- signif(assayDataElement(object, "probdloss"), digits=3)
      colnames(probdel) <- paste0("probdel.", colnames(probdel))
      df <- cbind(df, probdel)
    }
    probloss <- signif(assayDataElement(object, "probloss"), digits=3)
    colnames(probloss) <- paste0("probloss.", colnames(probloss))
    df <- cbind(df, probloss)
    probnorm <- signif(assayDataElement(object, "probnorm"), digits=3)
    colnames(probnorm) <- paste0("probnorm.", colnames(probnorm))
    df <- cbind(df, probnorm)
    probgain <- signif(assayDataElement(object, "probgain"), digits=3)
    colnames(probgain) <- paste0("probgain.", colnames(probgain))
    df <- cbind(df, probgain)
    if ("probamp" %in% assayDataElementNames(object)) {
      probamp <- signif(assayDataElement(object, "probamp"), digits=3)
      colnames(probamp) <- paste0("probamp.", colnames(probamp))
      df <- cbind(df, probamp)
    }
  }
  df
}

setMethod("plot", signature(x="cghRaw", y="missing"),
  getMethod("plot", signature=c(x="QDNAseqSignals", y="missing")))
setMethod("plot", signature(x="cghSeg", y="missing"),
  getMethod("plot", signature=c(x="QDNAseqSignals", y="missing")))
setMethod("plot", signature(x="cghCall", y="missing"),
  getMethod("plot", signature=c(x="QDNAseqSignals", y="missing")))

setMethod("frequencyPlot", signature(x="cghCall", y="missing"),
  getMethod("frequencyPlot", signature=c(x="QDNAseqCopyNumbers", y="missing")))

# EOF
