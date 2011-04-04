# ANALYSIS Visualisation/"Idiogram" (Plots the idiogram of genes on all chromosomes in the selected list. 
# For each chromosome, the Y axis reports the position of the gene along the cytobands. On the X axis the 
# fold changes of each gene are reported. Up-regulated genes are highlighted with red color, the 
# down-regulated with green color. Currently, this works only for human, mouse and rat data.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT idiogram.pdf
# PARAMETER chip.to.plot INTEGER FROM 1 TO 1000 DEFAULT 1 (Which chip to plot)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Plotting the idiogram of expressed genes
# DG 16.10.2006
# Modified by JTT 22.11.2006

# Loads the libraries
library(idiogram)
library(geneplotter)
library(genefilter)

# Loads the human idiogram data
data(Hs.cytoband)

# Parameter settings (default) for testing purposes
#chip.to.plot<-1
#image.width<-600
#image.height<-600

# Renaming variables
w<-image.width
h<-image.height
chip<-chip.to.plot

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Scaling the data to the same mean
scaled.dat<-genescale(dat2)

# Creating locations of genes
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
if(phenodata$chiptype[1]!="cDNA" | phenodata$chiptype[1]!="Illumina") {
   # Saves the chiptype into object lib
   lib<-phenodata$chiptype[1]
   lib<-as.character(lib)
}

# Account for the fact that annotation packages are from version 2.3 of Bioconductor
# named with an ".db" suffix. Add the suffix when missing to support data files
# from Chipster 1.3 and earlier. 
if (length(grep(".db", lib)) == 0 & length(grep("pmcdf", lib)) == 0) {
        lib <- paste(lib, ".db", sep="")
}

# Needs a parameter lib that defines the Affymetrix chip type
chromloc<-buildChromLocation(lib)

# Fold changes are stored in a named vector
fc<-scaled.dat[,chip]
names(fc)<-rownames(scaled.dat)

# Creating colors
cols<-rep("black", times=length(fc))
cols[fc>0]<-"red"
cols[fc<0]<-"green"

# What chromosomes are there for the species?
x<-get(paste(gsub(".db", "", lib), "CHR", sep=""))
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
chr<-unique(unlist(xx))
chr<-chr[-which(chr=="Un")]

# Fixing the chromosome locations object
chromloc@chromLocs<-chromloc@chromLocs[names(chromloc@chromLocs) %in% chr]
chromloc@chromInfo<-chromloc@chromInfo[names(chromloc@chromInfo) %in% chr]

# Plots the idiogram
pdf(file="idiogram.pdf", width=w/72, height=h/72)
midiogram(fc, chromloc, col=cols, pch=20) 
dev.off()