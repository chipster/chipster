# TOOL acgh-call-aberrations.R: "Call aberrations from segmented copy number data" (Call copy number aberrations from aCGH log ratios or NGS data.)
# INPUT segmented.tsv: segmented.tsv TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT aberrations.tsv: aberrations.tsv
# OUTPUT aberration-frequencies.pdf: aberration-frequencies.pdf
# PARAMETER number.of.copy.number.states: "Number of copy number states" TYPE [3: 3, 4: 4] DEFAULT 3 (Three states means calling loss vs. normal vs. gain, four states calls amplifications separately, and five also homozygous deletions.)
# PARAMETER genome.build: "Genome build" TYPE [GRCh37: GRCh37, NCBI36: NCBI36, NCBI35: NCBI35, NCBI34: NCBI34, other: other] DEFAULT GRCh37 (The genome build to use. GRCh37 = hg19, NCBI36 = hg18, NCBI35 = hg17, NCBI34 = hg16.)
# PARAMETER column: Cellularity TYPE METACOLUMN_SEL DEFAULT EMPTY (If available, phenodata column containing the cellularities (tumor cell percentages\) of the samples. If there are values larger than 1, they are assumed to be percentages and are divided by 100. Missing values are replaced by the mean value.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-06-27

source(file.path(chipster.common.path, "library-Chipster.R"))
source(file.path(chipster.common.path, "library-CGHcall.R"))

input <- readData("segmented.tsv")
phenodata <- readPhenodata("phenodata.tsv")
cgh <- toCgh(input, level="segmented")

organism <- "human"
if (genome.build == "other")
  organism <- "other"

cellularity <- rep(1, ncol(cgh))
if (column %in% colnames(phenodata))
  cellularity <- phenodata[,column]
if (max(cellularity, na.rm=TRUE) > 1)
  cellularity <- cellularity / 100
cellularity[is.na(cellularity)] <- mean(cellularity, na.rm=TRUE)
if (all(is.na(cellularity)))
  cellularity <- rep(1, ncol(cgh))

cgh.cal <- CGHcall(cgh, nclass=as.integer(number.of.copy.number.states), organism=organism, cellularity=cellularity, build=genome.build, ncpus=4)
cgh <- ExpandCGHcall(cgh.cal, cgh)

output <- fromCgh(cgh)
output <- addAnnotationColumns(input, output)
writeData(output, "aberrations.tsv")

# bitmap(file="aberration-frequencies.png", width=image.width/72, height=image.height/72)
pdf(file="aberration-frequencies.pdf", paper="a4r", width=0, height=0)
frequencyPlot(cgh)
dev.off()

# EOF
