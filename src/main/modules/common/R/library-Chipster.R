# Collection of handy code snippets to be shared between tools.
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-01-17

### Basic I/O for data and phenodata files.

readData <- function(file, header=TRUE, sep="\t", quote="", row.names=1, as.is=TRUE, check.names=FALSE, comment.char="", ...) {
  read.table(file, header=header, sep=sep, quote=quote, row.names=row.names, as.is=as.is, check.names=check.names, comment.char=comment.char, ...)
}

readPhenodata <- function(file, header=TRUE, sep="\t", quote="", as.is=TRUE, check.names=FALSE, comment.char="", ...) {
  read.table("phenodata.tsv", header=header, sep=sep, as.is=as.is, check.names=check.names, comment.char=comment.char, ...)
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

addAnnotationColumns <- function(input, output) {
  inputAnnotationColumns <- annotationColumns(input)
  outputAnnotationColumns <- annotationColumns(output)
  data.frame(output[, outputAnnotationColumns, drop=FALSE], input[rownames(output), setdiff(inputAnnotationColumns, outputAnnotationColumns), drop=FALSE], output[, setdiff(colnames(output), outputAnnotationColumns), drop=FALSE])
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
  chromosome[chromosome == "X"] <- "23"
  chromosome[chromosome == "Y"] <- "24"
  chromosome[chromosome == "MT"] <- "25"
  as.integer(chromosome)
}

# EOF
