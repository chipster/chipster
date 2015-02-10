# TOOL cna-segment.R: "Segment copy number data" (Segment copy number data before calling aberrations.)
# INPUT normalized-counts.tsv: normalized-counts.tsv TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT segmented-counts.tsv: segmented-counts.tsv
# PARAMETER minimum.number.of.probes.per.segment: "Minimum number of features per segment" TYPE [2: 2, 3: 3, 4: 4, 5: 5] DEFAULT 2 (Minimum number of features per segment.)
# PARAMETER minimum.number.of.sds.between.segments: "Minimum number of standard deviations between segments" TYPE DECIMAL FROM 0 TO 10 DEFAULT 1.0 (Minimum number of standard deviations required between segments.)
# PARAMETER significance.level: "Significance level" TYPE DECIMAL FROM 0 TO 1 DEFAULT 1E-10 (Significance level for the test to accept changepoints.)
# PARAMETER smoothBy: "Number of bins to smooth over" TYPE INTEGER FROM 1 TO 1000 DEFAULT 1 (For noisy data, the signal can be smoothed over n bins to reduce over-segmentation.)
# PARAMETER reNormalize: "Re-normalize" TYPE [TRUE: yes, FALSE: no] DEFAULT TRUE (Whether to perform an additional normalization step, which recursively searches for the interval containing the most segmented data, decreasing the interval length in each recursion. The recursive search makes the post-segmentation normalization robust against local maxima. This function is particularly useful for profiles for which, after segmentation, the 0-level does not coincide with many segments. It is more or less harmless to other profiles.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-23

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-QDNAseq.R'))

input <- readData("normalized-counts.tsv")
phenodata <- readPhenodata("phenodata.tsv")
copyNumbers <- toQDNAseqCopyNumbers(input, chiptype=phenodata$chiptype, level="copynumber")

copyNumbers <- segmentBins(copyNumbers, smoothBy=smoothBy, alpha=significance.level, min.width=as.integer(minimum.number.of.probes.per.segment), undo.SD=minimum.number.of.sds.between.segments)

if (reNormalize == "TRUE")
  copyNumbers <- normalizeSegmentedBins(copyNumbers)

output <- fromQDNAseqCopyNumbers(copyNumbers)
output <- addAnnotationColumns(input, output)
writeData(output, "segmented-counts.tsv")

# EOF
