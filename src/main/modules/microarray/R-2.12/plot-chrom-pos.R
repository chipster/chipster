# TOOL plot-chrom-pos.R: "Chromosomal position" (Plots the chromosomal positions of genes in the selected list. Currently, this works only for human, mouse and rat data. Before plotting, chips are scaled in order to infer the up- or down-regulation status.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT chromloc.pdf: chromloc.pdf 
# PARAMETER chip.to.plot: chip.to.plot TYPE INTEGER FROM 1 TO 1000 DEFAULT 1 (Which chip to plot)
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Plotting the chromosomal position of expressed genes
# JTT 13.6.2006
# Modified by JTT on 22.11.2006 according to DG

# Parameter settings (default) for testing purposes
#chip.to.plot<-c("chip.microarray001.cel")
#image.width<-c(600)
#image.height<-c(600)

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

# Account for the fact that annotation packages are from version 2.3 of Bioconductor
# named with an ".db" suffix. Add the suffix when missing to support data files
# from Chipster 1.3 and earlier. 
if (length(grep(".db", lib)) == 0 & length(grep("pmcdf", lib)) == 0) {
        lib <- paste(lib, ".db", sep="")
}

# Need a parameter lib that defines the Affymetrix chip type
chromloc<-buildChromLocation(gsub(".db", "", lib))

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Scaling the data to the same mean
scaled.dat<-genescale(dat2)

# Which genes are up- or down-regulated
up<-dat2[which(scaled.dat[,chip]>0),]
down<-dat2[which(scaled.dat[,chip]<0),]

# What chromosomes are there for the species?
x<-get(paste(gsub(".db", "", lib), "CHR", sep=""))
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
chr<-unique(unlist(xx))

# Plotting the chromosomes and genes
pdf(file="chromloc.pdf", width=w/72, height=h/72)
cPlot(chromloc, bg="White", fg="LightGrey", useChroms=chr)
cColor(rownames(up), "red", chromloc)
cColor(rownames(down), "green", chromloc)
dev.off()
