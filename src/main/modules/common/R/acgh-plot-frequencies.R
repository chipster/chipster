# TOOL acgh-plot-frequencies.R: "Plot copy number aberration frequencies" (Plot copy number aberration frequencies for different groups.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT group-frequencies.pdf: group-frequencies.pdf 
# PARAMETER column: Column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column defining sample groups.)
# PARAMETER chromosomes: Chromosomes TYPE STRING DEFAULT 0 (The numbers of the chromosomes to be plotted, separated by commas. 0 means all chromosomes. Ranges are also supported (e.g. 1,3,7-10\).)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-23

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-CGHregions.R'))

input <- readData("aberrations.tsv")
phenodata <- readPhenodata("phenodata.tsv")

cgh <- toCghRegions(input)

chrs.to.plot <- parseChromosomesToPlot(chromosomes, chromosomes(cgh))

# plot
pdf(file='group-frequencies.pdf', paper='a4r', width=0, height=0)
if (column == 'EMPTY') {
  frequencyPlot(cgh[chrs.to.plot, ])
} else {
  for (group in sort(unique(phenodata[, column]))) {
    if (is.na(group)) {
      frequencyPlot(cgh[chrs.to.plot, is.na(phenodata[,column])], main=paste(column, ' = ', group, sep=''))
    } else {
      frequencyPlot(cgh[chrs.to.plot, !is.na(phenodata[,column]) & phenodata[,column] == group], main=paste(column, ' = ', group, sep=''))
    }
  }
}
dev.off()

# EOF
