# TOOL acgh-normalize.R: "Normalize copy number data" (Normalize copy number data before segmentation.)
# INPUT raw.tsv: raw.tsv TYPE GENE_EXPRS 
# OUTPUT normalized.tsv: normalized.tsv
# PARAMETER method: Method TYPE [median: median, mode: mode, none: none] DEFAULT median (Normalization method.)
# PARAMETER smoothOutliers: "Smooth outliers" TYPE [TRUE: yes, FALSE: no] DEFAULT TRUE (Whether to smooth outliers, which is recommended before segmentation.)
# PARAMETER number.of.chromosomes: "Number of chromosomes" TYPE INTEGER DEFAULT 23 (Number of chromosomes. Usually 23 for sex-matched reference samples and 22 otherwise.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2013-06-27

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-CGHcall.R'))

input <- readData("raw.tsv")
cgh <- toCgh(input, level="copynumber")

cgh <- preprocess(cgh, nchrom=number.of.chromosomes)
cgh <- normalize(cgh, method=method, smoothOutliers=as.logical(smoothOutliers))

output <- fromCgh(cgh)
output <- addAnnotationColumns(input, output)
writeData(output, "normalized.tsv")

# EOF
