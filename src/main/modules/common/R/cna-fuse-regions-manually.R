# TOOL cna-fuse-regions-manually.R: "Fuse regions manually" (Fuses together neighboring regions that overlap an user-specified area of the genome.)
# INPUT regions.tsv: regions.tsv TYPE GENERIC 
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# OUTPUT fused-regions.tsv: regions.tsv 
# OUTPUT fused-region-frequencies.pdf: region-frequencies.pdf 
# PARAMETER position: Position TYPE STRING (Chromosomal region to fuse. Must contain three values that are separated by tabs, hyphens, colons or two dots (e.g. X:100-200 or 7:600..700\).)

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

# parse coordinates
position <- gsub('(\t|:|-|\\.\\.)', ';', position)
items <- strsplit(position, ';')[[1]]
chromosome <- items[1]
start <- as.integer(items[2])
end <- as.integer(items[3])

# match regions overlapping the area in question
condition <- input1$chromosome == chromosome & input1$start <= end & input1$end >= start

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
