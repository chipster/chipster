# ANALYSIS Visualisation/"Chromosomal position" (Plots the chromosomal positions of genes in the selected list. 
# Currently, this works only for human, mouse and rat data. Before plotting, chips are scaled in order to
# infer the up- or down-regulation status.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT chromloc.png
# PARAMETER chip.to.plot INTEGER FROM 1 TO 1000 DEFAULT 1 (Which chip to plot)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Plotting the chromosomal position of expressed genes
# JTT 13.6.2006
# Modified by JTT on 22.11.2006 according to DG

# Loads the libraries
library(geneplotter)
library(genefilter)

# Renaming variables
chip<-chip.to.plot
w<-image.width
h<-image.height

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Creating locations of genes
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
if(phenodata$chiptype[1]!="cDNA" | phenodata$chiptype[1]!="Illumina") {
   # Saves the chiptype into object lib
   lib<-phenodata$chiptype[1]
   lib<-as.character(lib)
}

# Need a parameter lib that defines the Affymetrix chip type
chromloc<-buildChromLocation(lib)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Scaling the data to the same mean
scaled.dat<-genescale(dat2)

# Which genes are up- or down-regulated
up<-dat2[which(scaled.dat[,chip]>0),]
down<-dat2[which(scaled.dat[,chip]<0),]

# Plotting the chromosomes and genes
bitmap(file="chromloc.png", width=w/72, height=h/72)
cPlot(chromloc, bg="White", fg="LightGrey")
cColor(rownames(up), "red", chromloc)
cColor(rownames(down), "green", chromloc)
dev.off()
