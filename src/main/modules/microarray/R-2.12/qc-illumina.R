# TOOL qc-illumina.R: Illumina (Illumina quality control using boxplots and density plots. This tool should be run on normalized data.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT boxplot.png: boxplot.png 
# OUTPUT densityplot.pdf: densityplot.pdf 
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Quality control of Illumina chips
# 30.7.2007 JTT
# 19.4.2011 MG, added legend to histogram

# Loads the libraries
library(limma)

# Renaming variables
w<-image.width
h<-image.height

# Reads in the data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Producing some basic plots of the data

# Boxplot
bitmap(file="boxplot.png", width=w/72, height=h/72)
boxplot(as.data.frame(dat2), las=2, cex.axis=0.5)
dev.off()

# Density plot
pdf(file="densityplot.pdf", width=w/72, height=h/72)
x<-c()
y<-c()
for(i in 1:ncol(dat2)) {
   x<-cbind(x, density(dat2[,i])$x)
   y<-cbind(y, density(dat2[,i])$y)
}
plot(x[,1], y[,1], col=1, ylim=c(0, max(y)), main="Intensities", type="l", lwd=2, lty=1, xlab = "log2(I)", ylab = "Density")
for(i in 2:ncol(dat2)) {
   lines(x[,i], y[,i], col=i, lwd=2, lty=1)
}
legend(legend=colnames(dat2), x="topright", lty=1, cex=0.75, col=c(1:ncol(dat2)))
dev.off()
