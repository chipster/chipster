# TOOL cna-normalize.R: "Normalize copy number data" (Normalizes copy number data before segmentation.)
# INPUT corrected-counts.tsv: corrected-counts.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT normalized-counts.tsv: normalized.tsv
# PARAMETER method: Method TYPE [median: median, mode: mode, none: none] DEFAULT median (Normalization method.)
# PARAMETER smoothOutliers: "Smooth outliers" TYPE [TRUE: yes, FALSE: no] DEFAULT TRUE (Whether to smooth outliers, which is recommended before segmentation.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-22

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-QDNAseq.R'))

input <- readData("corrected-counts.tsv")
phenodata <- readPhenodata("phenodata.tsv")
copyNumbers <- toQDNAseqCopyNumbers(input, chiptype=phenodata$chiptype, level="copynumber")

copyNumbers <- normalizeBins(copyNumbers, method=method)

if (smoothOutliers == "TRUE")
  copyNumbers <- smoothOutlierBins(copyNumbers)

output <- fromQDNAseqCopyNumbers(copyNumbers)
output <- addAnnotationColumns(input, output)
writeData(output, "normalized-counts.tsv")

# EOF
