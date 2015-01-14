# TOOL acgh-identify-regions.R: "Identify common regions from called copy number data" (Reduces dimensionality of called copy number data by identifying common breakpoints.)
# INPUT aberrations.tsv: aberrations.tsv TYPE GENE_EXPRS 
# OUTPUT regions.tsv: regions.tsv 
# OUTPUT region-frequencies.pdf: region-frequencies.pdf 
# PARAMETER max.info.loss: "Maximum loss of information" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.01 (Maximal information loss allowed.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-23

source(file.path(chipster.common.path, 'library-Chipster.R'))
source(file.path(chipster.common.path, 'library-CGHcall.R'))
source(file.path(chipster.common.path, 'library-CGHregions.R'))

input <- readData("aberrations.tsv")
cgh <- toCgh(input, level="calls")

regions <- CGHregions(cgh, max.info.loss)
regions2 <- regioning(cgh, max.info.loss, regions)

output <- fromRegioning(regions2)
writeData(output, "regions.tsv")

pdf(file="region-frequencies.pdf", paper="a4r", width=0, height=0)
frequencyPlot(regions)
plot(regions)
dev.off()

# EOF
