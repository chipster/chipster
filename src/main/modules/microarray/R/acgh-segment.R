# TOOL acgh-segment.R: "Segment copy number data" (Segment copy number data before calling aberrations.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT segmented.tsv: segmented.tsv 
# PARAMETER minimum.number.of.probes.per.segment: "Minimum number of features per segment" TYPE [2: 2, 3: 3, 4: 4, 5: 5] DEFAULT 2 (Minimum number of features per segment.)
# PARAMETER minimum.number.of.sds.between.segments: "Minimum number of standard deviations between segments" TYPE DECIMAL FROM 0 TO 10 DEFAULT 0 (Minimum number of standard deviations required between segments.)
# PARAMETER significance.level: "Significance level" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.01 (Significance level for the test to accept changepoints.)
# PARAMETER reNormalize: "Re-normalize" TYPE [TRUE: yes, FALSE: no] DEFAULT TRUE (Whether to perform an additional normalization step, which recursively searches for the interval containing the most segmented data, decreasing the interval length in each recursion. The recursive search makes the post-segmentation normalization robust against local maxima. This function is particularly useful for profiles for which, after segmentation, the 0-level does not coincide with many segments. It is more or less harmless to other profiles.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-06-27

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-CGHcall.R'))

input <- readData("normalized.tsv")
cgh <- toCgh(input, level="copynumber")

if (any(is.na(copynumber(cgh))))
  stop("CHIPSTER-NOTE: Data contains missing values. Please first run the tool Normalize copy number data.")

cgh <- segmentData(cgh, alpha=significance.level, min.width=as.integer(minimum.number.of.probes.per.segment), undo.SD=minimum.number.of.sds.between.segments)

if (reNormalize == "TRUE")
  cgh <- postsegnormalize(cgh)

output <- fromCgh(cgh)
output <- addAnnotationColumns(input, output)
writeData(output, "segmented.tsv")

# EOF
