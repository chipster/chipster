# TOOL acgh-cluster.R: "Cluster called copy number data" (Perform clustering of copy number data.)
# INPUT regions.tsv: regions.tsv TYPE GENERIC 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT wecca.pdf: wecca.pdf 
# PARAMETER type.of.calls: "Cluster on" TYPE [hard: "hard calls", soft: "call probabilities / soft calls", none: "none"] DEFAULT soft (Whether to cluster the arrays based on soft or hard calls. Hard calls are losses, normals, and gains, whereas soft calls refer to the respective probabilities of these calls. The preferred choice is to use soft calls whenever they are available. Choosing "none" will leave samples unordered.)
# PARAMETER column: Column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column to include in the output plot.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-01-25

source(file.path(chipster.common.path, "library-Chipster.R"))
source(file.path(chipster.common.path, "library-CGHregions.R"))

input <- readData("regions.tsv")
phenodata <- readData("phenodata.tsv")
regions <- toRegioning(input)

margins <- c(10, 1)
if (type.of.calls == "hard") {
  dendrogram <- WECCAhc(regions)
} else if (type.of.calls == "soft") {
  dendrogram <- WECCAsc(regions)
} else {
  dendrogram <- NA
  margins <- c(10, 2)
}

pdf(file="wecca.pdf", paper="a4", width=0, height=0)
if (column == "EMPTY") {
  WECCA.heatmap(regions, dendrogram, margins=margins)
} else {
  WECCA.heatmap(regions, dendrogram, margins=margins, ColSideColors=palette()[as.factor(phenodata[,column])])
}
dev.off()

# EOF
