# Construction end deconstruction of ExpressionSet objects.
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-23

library(Biobase)

toExpressionSet <- function(df, requirePositions=FALSE) {
  if (requirePositions) {
    pos <- c("chromosome", "start", "end")
    if (length(setdiff(pos, colnames(df))) != 0)
      stop("CHIPSTER-NOTE: The input file is required to contain the following columns: chromosome, start, end.")
    df$chromosome <- chromosomeToInteger(df$chromosome)
    df <- df[!is.na(df$chromosome) & !is.na(df$start) & !is.na(df$end), ]
    df <- df[order(df$chromosome, df$start, df$end), ]
    fData <- df[, pos]
    fData <- cbind(fData, df[, setdiff(annotationColumns(df), pos), drop=FALSE])
    colnames(fData)[1:3] <- c("Chromosome", "Start", "End")
  } else {
    fData <- df[, annotationColumns(df), drop=FALSE]
  }

  exprs <- as.matrix(df[, grep("^chip\\.", colnames(df)), drop=FALSE])
  samples <- sub("^chip\\.", "", colnames(exprs))
  colnames(exprs) <- samples
  new("ExpressionSet", exprs=exprs, featureData=new("AnnotatedDataFrame", data=fData))
}

# EOF
