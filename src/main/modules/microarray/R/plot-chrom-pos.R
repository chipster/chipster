# TOOL plot-chrom-pos.R: "Chromosomal position" (Plots the chromosomal positions of genes in the selected list. Currently, this works only for human, mouse and rat data. Before plotting, chips are scaled in order to infer the up- or down-regulation status.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT chromloc.pdf: chromloc.pdf 
# PARAMETER chip.to.plot: "Column to plot" TYPE COLUMN_SEL DEFAULT EMPTY (Column that contains the expression values. Data is mean centered)
# PARAMETER OPTIONAL image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER OPTIONAL image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# PARAMETER chip.to.plot: chip.to.plot TYPE INTEGER FROM 1 TO 1000 DEFAULT 1 (Which columns to plot)

# JTT 13.06.2006, created
# JTT 22.11.2006, modified according to DG
# MK 01.10.2013, small polishing
# MK 22.02.2014, added ability select column for plotting

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
#dat2<-dat[,grep("chip", names(dat))]
#
# Scaling the data to the same mean
#scaled.dat<-genescale(dat2)
#
# Which genes are up- or down-regulated
#if(chip > ncol(scaled.dat)) {
#	stop("CHIPSTER-NOTE: You have selected a chip that does not exists")
#}
#up<-dat2[which(scaled.dat[,chip]>0),]
#down<-dat2[which(scaled.dat[,chip]<0),]

if(is.numeric(dat[, grep(chip, colnames(dat))]) == FALSE) {
	stop("CHIPSTER-NOTE: You have selected a column that has no numeric entries")
}

#scaled.dat <-dat[,grep(chip, colnames(dat))]
if(length(grep("chip", chip)) > 0) {
	dat2 <- dat[,grep("chip", names(dat))]
	chip <- grep(chip, colnames(dat2))
	scaled.dat <- genescale(dat2)
	up <- dat2[which(scaled.dat[, chip] > 0),]
	down <- dat2[which(scaled.dat[, chip] < 0),]
} else {
	scaled.dat <- dat[, grep(chip, colnames(dat))]	
	up <- dat[which(scaled.dat > 0), ]
	down <- dat[which(scaled.dat < 0), ]
}

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
