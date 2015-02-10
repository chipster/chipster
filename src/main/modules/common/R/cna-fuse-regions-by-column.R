# TOOL cna-fuse-regions-by-column.R: "Fuse regions by column value" (Fuses together neighboring regions with e.g. a significant p-value for comparisons between groups.)
# INPUT regions.tsv: regions.tsv TYPE GENERIC 
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# OUTPUT fused-regions.tsv: regions.tsv 
# OUTPUT fused-region-frequencies.pdf: region-frequencies.pdf 
# PARAMETER column: Column TYPE COLUMN_SEL DEFAULT pvalue (Data column to use)
# PARAMETER criteria: Is TYPE [smaller: "smaller than", smaller-equal: "smaller than or equal to", equal: "equal to", larger-equal: "larger than or equal to", larger: "larger than", not-equal: "not equal to", within: within, outside: outside] DEFAULT smaller (Whether to match regions where column values that are smaller, larger, or equal compared to the cutoff. Use within and outside for symmetrical comparisons around zero.)
# PARAMETER cutoff: "Cut-off" TYPE DECIMAL DEFAULT 0.05 (Cutoff for matching.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-06-27

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-CGHcall.R'))
source(file.path(chipster.common.path, 'library-CGHregions.R'))

input1 <- readData("regions.tsv")
regions <- toCghRegions(input1)

input2 <- readData("aberrations.tsv")
cgh <- toCgh(input2, level="calls")

if (!identical(sampleNames(regions), sampleNames(cgh))) {
  if (length(setdiff(sampleNames(regions), sampleNames(cgh))) == 0) {
    cgh <- cgh[,sampleNames(regions)]
  } else {
    stop("CHIPSTER-NOTE: Sample names do not match.\nregions.tsv: ", paste(sampleNames(regions), collapse=", "), "\naberrations.tsv: ", paste(sampleNames(cgh), collapse=", "))
  }
}

if (criteria == 'smaller') {
  condition <- input1[,column] < cutoff
} else if (criteria == 'smaller-equal') {
  condition <- input1[,column] <= cutoff
} else if (criteria == 'equal') {
  condition <- input1[,column] == cutoff
} else if (criteria == 'larger-equal') {
  condition <- input1[,column] >= cutoff
} else if (criteria == 'larger') {
  condition <- input1[,column] > cutoff
} else if (criteria == 'not-equal') {
  condition <- input1[,column] != cutoff
} else if (criteria == 'within') {
  condition <- abs(input1[,column]) <= abs(cutoff)
} else if (criteria == 'outside') {
  condition <- abs(input1[,column]) > abs(cutoff)
}

prev.condition <- c(FALSE, condition[-length(condition)])
next.condition <- c(condition[-1], FALSE)

# new regions start where either:
# 1. condition does not match
# 2. previous condition did not match
# 3. previous region was in a different chromosome
starts <- !condition | !prev.condition |
  input1$chromosome != c(0, input1$chromosome[-nrow(input1)])

# new regions end where either:
# 1. condition does not match
# 2. next condition does not match
# 3. next region is in a different chromosome
ends <- !condition | !next.condition |
  input1$chromosome != c(input1$chromosome[-1], 0)

newreg <- input1[starts, c('chromosome', 'start', 'end')]
newreg$end <- input1$end[ends]
newRegions <- CGHregionsManual(cgh, newreg)
newRegioning <- regioning(cgh, NA, newRegions)

output <- fromRegioning(newRegioning)
writeData(output, "fused-regions.tsv")

pdf(file="fused-region-frequencies.pdf", paper="a4r", width=0, height=0)
frequencyPlot(newRegions)
plot(newRegions)
dev.off()

# EOF
