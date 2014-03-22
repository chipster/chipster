# TOOL cna-call-aberrations.R: "Call aberrations from segmented copy number data" (Call copy number aberrations from aCGH log ratios or NGS data.)
# INPUT segmented-counts.tsv: segmented-counts.tsv TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT aberrations.tsv: aberrations.tsv
# OUTPUT aberration-frequencies.pdf: aberration-frequencies.pdf
# PARAMETER number.of.copy.number.states: "number of copy number states" TYPE [3: 3, 4: 4] DEFAULT 3 (Three states means calling loss vs. normal vs. gain, four states calls amplifications separately, and five also homozygous deletions.)
# PARAMETER column: cellularity TYPE METACOLUMN_SEL DEFAULT EMPTY (If available, phenodata column containing the cellularities (tumor cell percentages\) of the samples. If there are values larger than 1, they are assumed to be percentages and are divided by 100. Missing values are replaced by the mean value.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-23

source(file.path(chipster.common.path, "library-Chipster.R"))
source(file.path(chipster.common.path, "library-QDNAseq.R"))

input <- readData("segmented-counts.tsv")
phenodata <- readPhenodata("phenodata.tsv")
copyNumbers <- toQDNAseqCopyNumbers(input, chiptype=phenodata$chiptype, level="segmented")

organism <- "human"
build <- sub("\\..*", "", phenodata$chiptype[1])
if (substr(build, 1, 2) != "hg")
  organism <- "other"

cellularity <- rep(1, ncol(copyNumbers))
if (column %in% colnames(phenodata))
  cellularity <- phenodata[,column]
if (max(cellularity, na.rm=TRUE) > 1)
  cellularity <- cellularity / 100
cellularity[is.na(cellularity)] <- mean(cellularity, na.rm=TRUE)
if (all(is.na(cellularity)))
  cellularity <- rep(1, ncol(copyNumbers))

copyNumbers <- callBins(copyNumbers, nclass=as.integer(number.of.copy.number.states), organism=organism, cellularity=cellularity, build=build, ncpus=4)

output <- fromQDNAseqCopyNumbers(copyNumbers)
output <- addAnnotationColumns(input, output)
writeData(output, "aberrations.tsv")

# bitmap(file="aberration-frequencies.png", width=image.width/72, height=image.height/72)
pdf(file="aberration-frequencies.pdf", paper="a4r", width=0, height=0)
frequencyPlot(copyNumbers, main=paste0("Frequency Plot (n=", ncol(copyNumbers), ")"))
dev.off()

# EOF
