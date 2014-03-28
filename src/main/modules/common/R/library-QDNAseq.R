# Construction end deconstruction of QDNAseq objects.
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-23

library(Biobase)
library(QDNAseq)

toQDNAseqReadCounts <- function(df, chiptype) {
  if (length(unique(chiptype)) != 1)
    stop("CHIPSTER-NOTE: Umabiguous chiptype definition:\n", paste(unique(chiptype), collapse=", "))
  if (!file.exists(file.path(chipster.tools.path, "QDNAseq", paste0("QDNAseq.", chiptype[1], ".rds"))))
    stop("CHIPSTER-NOTE: Unknown chiptype:\n", chiptype[1])
  bins <- readRDS(file.path(chipster.tools.path, "QDNAseq", paste0("QDNAseq.", chiptype[1], ".rds")))
  chip <- as.matrix(df[, grep("^chip\\.", colnames(df)), drop=FALSE])
  samples <- sub("^chip\\.", "", colnames(chip))
  m <- matrix(nrow=nrow(bins), ncol=ncol(chip), dimnames=list(featureNames(bins), samples))
  m[rownames(input),] <- chip
  bins@data$use <- FALSE
  bins@data[rownames(df), "use"] <- TRUE
  object <- new("QDNAseqReadCounts", bins=bins, counts=m, phenodata=data.frame(samples, row.names=samples))
  object
}

toQDNAseqCopyNumbers <- function(df, chiptype, level="any") {
  if (length(unique(chiptype)) != 1)
    stop("CHIPSTER-NOTE: Umabiguous chiptype definition:\n", paste(unique(chiptype), collapse=", "))
  if (!file.exists(file.path(chipster.tools.path, "QDNAseq", paste0("QDNAseq.", chiptype[1], ".rds"))))
    stop("CHIPSTER-NOTE: Unknown chiptype:\n", chiptype[1])
  bins <- readRDS(file.path(chipster.tools.path, "QDNAseq", paste0("QDNAseq.", chiptype[1], ".rds")))
  chip <- as.matrix(df[, grep("^chip\\.", colnames(df)), drop=FALSE])
  samples <- sub("^chip\\.", "", colnames(chip))
  m <- matrix(nrow=nrow(bins), ncol=ncol(chip), dimnames=list(featureNames(bins), samples))
  m[rownames(input),] <- chip
  bins@data$use <- FALSE
  bins@data[rownames(df), "use"] <- TRUE
  object <- new("QDNAseqCopyNumbers", bins=bins, copynumber=QDNAseq:::unlog2adhoc(m), phenodata=data.frame(samples, row.names=samples))
  if (level == "copynumber")
    return(object)
  segmented <- as.matrix(df[, grep("^segmented\\.", colnames(df)), drop=FALSE])
  if (ncol(segmented) == 0) {
    if (level == "any") {
      return(object)
    } else {
      stop("CHIPSTER-NOTE: The input file does not contain segmented data. Please segment it first.")
    }
  }
  m[rownames(input),] <- segmented
  assayDataElement(object, "segmented") <- QDNAseq:::unlog2adhoc(m)
  if (level == "segmented")
    return(object)
  calls <- as.matrix(df[, grep("^flag\\.", colnames(df)), drop=FALSE])
  if (ncol(calls) == 0) {
    if (level == "any") {
      return(object)
    } else {
      stop("CHIPSTER-NOTE: The input file does not contain called data. Please call aberrations first.")
    }
  }
  m[rownames(input),] <- calls
  assayDataElement(object, "calls") <- m
  probdel <- as.matrix(df[, grep("^probdel\\.", colnames(df)), drop=FALSE])
  if (ncol(probdel) > 0) {
    m[rownames(input),] <- probdel
    assayDataElement(object, "probdel") <- m
  }
  probloss <- as.matrix(df[, grep("^probloss\\.", colnames(df)), drop=FALSE])
  m[rownames(input),] <- probloss
  assayDataElement(object, "probloss") <- m
  probnorm <- as.matrix(df[, grep("^probnorm\\.", colnames(df)), drop=FALSE])
  m[rownames(input),] <- probnorm
  assayDataElement(object, "probnorm") <- m
  probgain <- as.matrix(df[, grep("^probgain\\.", colnames(df)), drop=FALSE])
  m[rownames(input),] <- probgain
  assayDataElement(object, "probgain") <- m
  probamp <- as.matrix(df[, grep("^probamp\\.", colnames(df)), drop=FALSE])
  if (ncol(probamp) > 0) {
    m[rownames(input),] <- probamp
    assayDataElement(object, "probamp") <- m
  }
  object
}

toQDNAseqSignals <- function(df, ...) {
  tmp <- as.matrix(df[, grep("^chip\\.", colnames(df)), drop=FALSE])
  if (storage.mode(tmp) == "integer") {
    signals <- toQDNAseqReadCounts(df, ...)
  } else {
    signals <- toQDNAseqCopyNumbers(df, ...)
  }
  signals
}

fromQDNAseqReadCounts <- function(object) {
  if ("use" %in% colnames(fData(object)))
    object <- object[fData(object)$use,]
  samples <- sampleNames(object)
  df <- fData(object)[,c("chromosome", "start", "end")]
  counts <- assayDataElement(object, "counts")
  colnames(counts) <- paste0("chip.", colnames(counts))
  df <- cbind(df, counts)
  df
}

fromQDNAseqCopyNumbers <- function(object) {
  if ("use" %in% colnames(fData(object)))
    object <- object[fData(object)$use,]
  samples <- sampleNames(object)
  df <- fData(object)[,c("chromosome", "start", "end")]
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
  copynumber <- signif(QDNAseq:::log2adhoc(assayDataElement(object, "copynumber")), digits=3)
  colnames(copynumber) <- paste0("chip.", colnames(copynumber))
  df <- cbind(df, copynumber)
  if ("segmented" %in% assayDataElementNames(object)) {
    segmented <- signif(QDNAseq:::log2adhoc(assayDataElement(object, "segmented")), digits=3)
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

# EOF
