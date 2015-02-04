# TOOL cna-correct.R: "Correct for GC content and mappability" (This tool corrects NGS data for GC content for copy number analysis.)
# INPUT filtered-counts.tsv: "Data table with read counts" TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT corrected-counts.tsv: "Data table with read counts corrected for GC content and mappability"
# PARAMETER method: Method TYPE [ratio: ratio, median: median, none: none] DEFAULT ratio (Whether to perform the correction as a ratio to the loess fit, or as an additive correction to the level of the median loess value.)
# PARAMETER span: "Degree of smoothing" TYPE DECIMAL DEFAULT 0.65 (The alpha/span parameter for the degree of smoothing, controls the neighboorhood used for the fitting. For alpha < 1, the neighbourhood includes proportion alpha of the points, and these have tricubic weighting. For alpha > 1, all points are used, with the ‘maximum distance’ assumed to be alpha^1/2.)
# PARAMETER family: Family TYPE [gaussian: gaussian, symmetric: symmetric] DEFAULT symmetric (Family used for the fitting. For gaussian, fitting is done by weighted least squares. For symmetric, a few iterations of an M-estimation procedure with Tukey's biweight are used.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-22

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-QDNAseq.R'))

input <- readData("filtered-counts.tsv")
phenodata <- readPhenodata("phenodata.tsv")
readCounts <- toQDNAseqReadCounts(input, chiptype=phenodata$chiptype)

copyNumbers <- correctBins(readCounts, method=method, span=span, family=family, force=TRUE)

output <- fromQDNAseqCopyNumbers(copyNumbers)
output <- addAnnotationColumns(input, output)
writeData(output, "corrected-counts.tsv")

# EOF
