# TOOL cna-filter-bins.R: "Filter copy number bins" (Filters out low-quality bins.)
# INPUT read-counts.tsv: read-counts.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT filtered-counts.tsv: filtered-counts.tsv
# PARAMETER qdnaseq: "QDNAseq blacklist" TYPE [TRUE: remove, FALSE: keep] DEFAULT TRUE (Whether to remove bins in the QDNAseq blacklist, which is based on 1000 Genomes data.)
# PARAMETER encode: "ENCODE blacklist" TYPE [TRUE: remove, FALSE: keep] DEFAULT TRUE (Whether to remove bins overlapping with the ENCODE blacklisted regions.)
# PARAMETER mappability: "Minimum mappability" TYPE INTEGER FROM 0 TO 100 DEFAULT 0 (Mappability threshold to filter out bins with mappabilities lower than the given value.)
# PARAMETER allosomes: "Sex chromosomes" TYPE [TRUE: remove, FALSE: keep] DEFAULT TRUE (Whether to filter out X and Y chromosomes.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-04-01

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-QDNAseq.R'))

input <- readData("read-counts.tsv")
phenodata <- readPhenodata("phenodata.tsv")
readCounts <- toQDNAseqReadCounts(input, chiptype=phenodata$chiptype)

readCounts <- applyFilters(readCounts, residual=as.logical(qdnaseq), blacklist=as.logical(encode), mappability=mappability, filterAllosomes=as.logical(allosomes))

output <- fromQDNAseqReadCounts(readCounts)
output <- addAnnotationColumns(input, output)
writeData(output, "filtered-counts.tsv")

# EOF
