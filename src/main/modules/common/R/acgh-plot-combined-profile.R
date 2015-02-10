# TOOL acgh-plot-combined-profile.R: "Plot profiles of matched copy number and expression" (Plot profiles of two priorly matched data sets of copy number and expression. This tool must be run on the output from the tool Match copy number and expression features.)
# INPUT matched-cn-and-expression.tsv: matched-cn-and-expression.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT matched-cn-and-expression-profile.pdf: matched-cn-and-expression-profile.pdf 
# PARAMETER samples: Samples TYPE STRING DEFAULT 1 (The numbers of the samples to be plotted, separated by commas. Ranges are also supported (e.g. 1,3,7-10\).)
# PARAMETER chromosome: Chromosome TYPE INTEGER DEFAULT 0 (The chromosome to plot. Use 0 for all.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-25

source(file.path(chipster.common.path, "library-Chipster.R"))
source(file.path(chipster.common.path, "library-intCNGEan.R"))

# read the input files
input <- readData("matched-cn-and-expression.tsv")
phenodata <- readPhenodata("phenodata.tsv")

matched <- toIntCNGEanMatched(input)
sampleNames(matched$CNdata.matched) <- phenodata$description
sampleNames(matched$GEdata.matched) <- phenodata$description

samples.to.plot <- parseSamplesToPlot(samples, 1:nrow(phenodata))

# plot
pdf(file="matched-cn-and-expression-profile.pdf", paper="a4r", width=0, height=0)
for (sample in samples.to.plot)
  intCNGEan.profilesPlot(matched$CNdata.matched, matched$GEdata.matched, sampleNo=sample, chr=chromosome)
dev.off()

# EOF
