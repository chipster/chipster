# TOOL cna-define-experiment.R: "Define CNA-seq experiment" (This tool counts reads in fixed-sized bins and creates a phenodata file containing descriptive information about samples and experiment setup.)
# INPUT alignment{...}.bam: "Input BAM files" TYPE GENERIC
# OUTPUT read-counts.tsv: "Data table with read counts"
# OUTPUT META phenodata.tsv: "Experiment description file"
# PARAMETER experiment: Experiment TYPE [hg19.1kbp.SR50: "genome=hg19, binSize=1kbp, type=SR50", hg19.5kbp.SR50: "genome=hg19, binSize=5kbp, type=SR50", hg19.10kbp.SR50: "genome=hg19, binSize=10kbp, type=SR50", hg19.15kbp.SR50: "genome=hg19, binSize=15kbp, type=SR50", hg19.30kbp.SR50: "genome=hg19, binSize=30kbp, type=SR50", hg19.50kbp.SR50: "genome=hg19, binSize=50kbp, type=SR50", hg19.100kbp.SR50: "genome=hg19, binSize=100kbp, type=SR50", hg19.500kbp.SR50: "genome=hg19, binSize=500kbp, type=SR50", hg19.1000kbp.SR50: "genome=hg19, binSize=1000kbp, type=SR50"] DEFAULT hg19.15kbp.SR50 (Experiment type, including organism, genome build, bin size to use, and the type of sequencing experiment (e.g. SR50 for single-read 50 bp\).)
# PARAMETER isFirstMateRead: "Is first mate read" TYPE [FALSE: no, TRUE: yes, NA: any] DEFAULT NA (Whether the first mate read should be returned or not, or whether mate read number should be ignored.)
# PARAMETER isSecondMateRead: "Is second mate read" TYPE [FALSE: no, TRUE: yes, NA: any] DEFAULT NA (Whether the second mate read should be returned or not, or whether mate read number should be ignored.)
# PARAMETER isNotPrimaryRead: "Is not primary read" TYPE [FALSE: no, TRUE: yes, NA: any] DEFAULT NA (Whether alignments that are primary, are not primary or whose primary status does not matter should be returned. A non-primary alignment (“secondary alignment” in the SAM specification\) might result when a read aligns to multiple locations. One alignment is designated as primary and has this flag set to FALSE; the remainder, for which this flag is TRUE, are designated by the aligner as secondary.)
# PARAMETER isNotPassingQualityControls: "Is not passing quality controls" TYPE [FALSE: no, TRUE: yes, NA: any] DEFAULT FALSE (Whether reads passing quality controls, reads not passing quality controls, or any read should be returned.)
# PARAMETER isDuplicate: "Is duplicate" TYPE [FALSE: no, TRUE: yes, NA: any] DEFAULT FALSE (Wether un-duplicated, duplicated, or any reads should be returned. 'Duplicated' reads may represent PCR or optical duplicates.)
# PARAMETER minMapq: "Minimum mapping quality" TYPE INTEGER DEFAULT 37 (Minimum mapping quality. Reads with a lower MAPQ are not counted.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-04-01

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-QDNAseq.R'))

bins <- readRDS(file.path(chipster.tools.path, 'QDNAseq', paste0('QDNAseq.', experiment, '.rds')))

readCounts <- binReadCounts(bins,
  isPaired=NA,
  isProperPair=NA,
  isUnmappedQuery=FALSE,
  hasUnmappedMate=NA,
  isMinusStrand=NA,
  isMateMinusStrand=NA,
  isFirstMateRead=as.logical(isFirstMateRead),
  isSecondMateRead=as.logical(isSecondMateRead),
  isNotPrimaryRead=as.logical(isNotPrimaryRead),
  isNotPassingQualityControls=as.logical(isNotPassingQualityControls),
  isDuplicate=as.logical(isDuplicate),
  minMapq=minMapq)

fData(readCounts)$use <- TRUE

# generate phenodata
phenodata <- data.frame(sample=paste0(readCounts$name, '.bam'), experiment='QDNAseq', chiptype=experiment, reads=readCounts$total.reads, group='', stringsAsFactors=FALSE)

# write outputs
output <- fromQDNAseqReadCounts(readCounts)
writeData(output, "read-counts.tsv")
writePhenodata(phenodata, "phenodata.tsv")

# EOF
