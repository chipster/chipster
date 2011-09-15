# ANALYSIS Visualisation/"Boxplot" (Creates a boxplot of normalized data. One box per chip is plotted.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv 
# OUTPUT boxplot.png
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Boxplot
# JTT 2.10.2007
# MG, 15.9.2011
# Updated coloring

# Parameter settings (default) for testing purposes
#image.width<-c(600)
#image.height<-c(600)

# Renaming variables
w<-image.width
h<-image.height

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Plotting
if(nrow(phenodata)==ncol(dat2)) {
   bitmap(file="boxplot.png", width=w/72, height=h/72)
   par(mar=c(12,5,5,5))
   boxplot(as.data.frame(dat2), las=2, names=phenodata$description, col=phenodata$group+2)
   dev.off()
} else {
   bitmap(file="boxplot.png", width=w/72, height=h/72)
   par(mar=c(12,5,5,5))
   boxplot(as.data.frame(dat2), las=2, col=phenodata$group+2)
   dev.off()
}

