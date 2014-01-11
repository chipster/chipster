# Collection of handy code snippets to be shared between tools.
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-01-11

# This function identifies matrices (chip, flag, segmented, ...) present in the
# data and returns the names of annotation columns, i.e. the ones that are not
# part of any matrix.
annotationColumns <- function(columns) {
  suffix <- sub("^chip\\.", "", columns[grep("^chip\\.", columns)[1]])
  suffix <- paste0(suffix, "$")
  matrices <- sub(suffix, "", columns[grep(suffix, columns)])
  annotations <- seq_along(columns)
  for (m in matrices)
    annotations <- setdiff(annotations, grep(m, columns))
  columns[annotations]
}

# EOF
