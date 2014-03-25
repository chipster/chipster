# Collection of handy code snippets to be shared between tools.
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-25

### Basic I/O for data and phenodata files.

readData <- function(file, header=TRUE, sep="\t", quote="", row.names=1, as.is=TRUE, check.names=FALSE, comment.char="", ...) {
  read.table(file, header=header, sep=sep, quote=quote, row.names=row.names, as.is=as.is, check.names=check.names, comment.char=comment.char, ...)
}

readPhenodata <- function(file, header=TRUE, sep="\t", quote="", as.is=TRUE, check.names=FALSE, comment.char="", ...) {
  read.table(file, header=header, sep=sep, as.is=as.is, check.names=check.names, comment.char=comment.char, ...)
}

writeData <- function(x, file, quote=FALSE, sep="\t", na="", ...) {
  options(scipen=10)
  write.table(x, file=file, quote=quote, sep=sep, na=na, ...)
}

writePhenodata <- function(x, file, quote=FALSE, sep="\t", na="", row.names=FALSE, ...) {
  options(scipen=10)
  write.table(x, file=file, quote=quote, sep=sep, na=na, row.names=row.names, ...)
}

### Miscallenous helper functions.

addAnnotationColumns <- function(input, output, rows=rownames(output), exclude=NULL) {
  inputAnnotationColumns <- annotationColumns(input)
  outputAnnotationColumns <- annotationColumns(output)
  outputDataColumns <- dataColumns(output)
  annotationColumnsToAdd <- setdiff(inputAnnotationColumns, outputAnnotationColumns)
  annotationColumnsToAdd <- setdiff(annotationColumnsToAdd, exclude)
  data.frame(output[, outputAnnotationColumns, drop=FALSE], input[rows, annotationColumnsToAdd, drop=FALSE], output[, outputDataColumns, drop=FALSE])
}

annotationColumns <- function(df) {
  columns <- colnames(df)
  suffix <- sub("^chip\\.", "", columns[grep("^chip\\.", columns)[1]])
  suffix <- paste0(suffix, "$")
  matrices <- sub(suffix, "", columns[grep(suffix, columns)])
  annotations <- seq_along(columns)
  for (m in matrices)
    annotations <- setdiff(annotations, grep(m, columns))
  columns[annotations]
}

chromosomeToCharacter <- function(chromosome) {
  chromosome <- as.character(chromosome)
  chromosome[chromosome == "23"] <- "X"
  chromosome[chromosome == "24"] <- "Y"
  chromosome[chromosome == "25"] <- "MT"
  chromosome
}

chromosomeToInteger <- function(chromosome) {
  chromosome <- toupper(chromosome)
  chromosome[chromosome == "X"] <- "23"
  chromosome[chromosome == "Y"] <- "24"
  chromosome[chromosome == "MT" | chromosome == "M"] <- "25"
  as.integer(chromosome)
}

dataColumns <- function(df) {
  columns <- colnames(df)
  suffix <- sub("^chip\\.", "", columns[grep("^chip\\.", columns)[1]])
  suffix <- paste0(suffix, "$")
  matrices <- sub(suffix, "", columns[grep(suffix, columns)])
  notAnnotations <- integer()
  for (m in matrices)
    notAnnotations <- c(notAnnotations, grep(m, columns))
  columns[notAnnotations]
}

parseChromosomesToPlot <- function(string, chromosomes) {
  if (class(chromosomes) == "character")
    chromosomes <- chromosomeToInteger(chromosomes)
  string <- gsub("X", "23", string, ignore.case=TRUE)
  string <- gsub("Y", "24", string, ignore.case=TRUE)
  string <- gsub("M[T]?", "25", string, ignore.case=TRUE)
  string <- gsub("[^0-9,-]", ",", string)
  items <- strsplit(string, ",")[[1]]
  chromosomesToPlot <- integer()
  for (item in items) {
    item <- item[item != ""]
    if (length(item) == 0) next
    range <- strsplit(item, "-")[[1]]
    range <- range[range != ""]
    if (length(range) == 0) next
    chromosomesToPlot <- c(chromosomesToPlot, seq(as.integer(range[1]), as.integer(range[length(range)])))
  }
  chromosomesToPlot <- unique(chromosomesToPlot)
  chromosomesToPlot <- chromosomesToPlot[chromosomesToPlot %in% chromosomes]
  if (length(chromosomesToPlot) == 0)
    chromosomesToPlot <- 0
  if (0 %in% chromosomesToPlot) {
    cond <- rep(TRUE, length(chromosomes))
  } else {
    cond <- chromosomes %in% chromosomesToPlot
  }
  cond
}

parseSamplesToPlot <- function(string, from) {
  if (string == "0")
    return(from)
  string <- gsub("[^0-9,-]", ",", string)
  items <- strsplit(string, ",")[[1]]
  samplesToPlot <- integer()
  for (item in items) {
    item <- item[item != ""]
    if (length(item) == 0) next
    range <- strsplit(item, "-")[[1]]
    range <- range[range != ""]
    if (length(range) == 0) next
    samplesToPlot <- c(samplesToPlot, seq(as.integer(range[1]), as.integer(range[length(range)])))
  }
  samplesToPlot <- unique(samplesToPlot)
  samplesToPlot <- samplesToPlot[samplesToPlot %in% from]
  if (length(samplesToPlot) == 0)
    return(from)
  samplesToPlot
}

# EOF
