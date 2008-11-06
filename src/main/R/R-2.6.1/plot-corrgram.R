# ANALYSIS Visualisation/"Correlogram" (Plots the correlations between samples in a graph.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT corrgram.png
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Correlogram
# JTT 18.10.2007

# Renaming variables
w<-image.width
h<-image.height

# Loads libraries
library(corrgram)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Plotting
bitmap(file="corrgram.png", width=w/72, height=h/72)
corrgram(x=dat2, lower.panel=panel.shade, upper.panel=panel.pts, pch=".")
dev.off()


