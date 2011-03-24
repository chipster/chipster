# ANALYSIS "aCGH"/"Cluster called aCGH data" (Perform clustering of aCGH arrays.)
# INPUT GENERIC regions.tsv, GENERIC phenodata.tsv
# OUTPUT wecca.pdf
# PARAMETER type.of.calls [hard, soft] DEFAULT soft (Whether to cluster the arrays based on soft or hard calls. Hard calls are losses, normals, and gains, whereas soft calls refer to the respective probabilities of these calls. The preferred choice is to use soft calls whenever they are available.)
# PARAMETER image.width INTEGER FROM 200 TO 6400 DEFAULT 2400 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 6400 DEFAULT 2400 (Height of the plotted network image)

# cluster-acgh.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2010-10-12

library(WECCA)

dat <- read.table('regions.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
phenodata <- read.table("phenodata.tsv", header=TRUE, sep="\t", as.is=TRUE)

dat$chromosome[dat$chromosome=='X'] <- '23'
dat$chromosome[dat$chromosome=='Y'] <- '24'
dat$chromosome[dat$chromosome=='MT'] <- '25'
dat$chromosome <- as.integer(dat$chromosome)

ann <- dat[,c('chromosome', 'start', 'end')]
colnames(ann) <- c('Chromosome', 'Start', 'End')
hardcalls <- as.matrix(dat[,grep("^flag", names(dat))])
if (ncol(hardcalls)==0)
  stop('CHIPSTER-NOTE: No copy number calls were found. Please run tools Call copy number aberrations from aCGH data and Identify common regions from called aCGH data first.')
colnames(hardcalls) <- phenodata$description
softcalls <- as.matrix(dat[,grep("^prob", names(dat))])
if (ncol(softcalls)==0) {
  if (type.of.calls == 'soft')
    stop('CHIPSTER-NOTE: No soft calls were found. Please try running with parameter type.of.calls=hard.')
  # the size of the softcalls matrix is used to detect whether calling was performed with 3 or 4 copy number states.
  # therefore we will construct such a matrix
  if (2 %in% hardcalls) {
    softcalls <- cbind(hardcalls, hardcalls, hardcalls, hardcalls)
  } else {
    softcalls <- cbind(hardcalls, hardcalls, hardcalls)
  }
}
colnames(softcalls) <- sub('\\.', '_', colnames(softcalls))

regions <- list(ann=ann, hardcalls=hardcalls, softcalls=softcalls)

if (type.of.calls == 'hard') {
  dendrogram <- WECCAhc(regions)
} else {
  dendrogram <- WECCAsc(regions)
}

# there's a small bug in the WECCA library, so here's a fixed version of this function:
WECCA.heatmap.modified <- function (cghdata.regioned, dendrogram) 
{
    nclass <- dim(cghdata.regioned$softcalls)[2]/dim(cghdata.regioned$hardcalls)[2]
    chr.color <- rep("blue", dim(cghdata.regioned$hardcalls)[1])
    ids <- ((cghdata.regioned$ann[, 1]%%2) == 0)
    chr.color[ids] <- c("yellow")
    Y <- rep(FALSE, dim(cghdata.regioned$hardcalls)[1])
    for (i in 2:(dim(cghdata.regioned$ann)[1])) {
        if ((cghdata.regioned$ann[i - 1, 1] != cghdata.regioned$ann[i, 1])) {
            Y[i] <- TRUE
        }
    }
    Y[1] <- TRUE
    begin.chr <- rep("", dim(cghdata.regioned$ann)[1])
    begin.chr[Y] <- cghdata.regioned$ann[Y, 1]
    color.coding <- c("red", "black", "green", "white")[1:nclass]
    heatmap(cghdata.regioned$hardcalls, Colv = as.dendrogram(dendrogram), 
        Rowv = NA, col = color.coding, labRow = begin.chr, RowSideColors = chr.color, 
        scale = "none")
}

pdf(file='wecca.pdf', width=image.width/72, height=image.height/72)
WECCA.heatmap.modified(regions, dendrogram)
dev.off()

# EOF
