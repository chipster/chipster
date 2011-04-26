# ANALYSIS Visualisation/"Volcano plot" (Tests whether the genes are differentially expressed using one-sample t-test,
# and plots the result in a form of a Volcano plot, but a table with the original data, the scaled fold change values and the adjusted p-values is also output.
# This tool is best suited for 2-color array data, where the comparison to a common reference is done on the array.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT volcanoP.pdf, volcanoSE.pdf, one-sample.tsv
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER expression.threshold DECIMAL FROM -10 TO 10 DEFAULT 1 (Expression threshold for plotting) 
# PARAMETER p.value.adjustment.method [none, Bonferroni, Holm, Hochberg, BH, BY] DEFAULT BH (Multiple testing correction method)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# Volcano plot
# JTT 19.9.2007
#
# MG. 10.2.2011
# Added table with FC and p-values to the output
# Added p-value adjustment capability

# Renaming variables
w<-image.width
h<-image.height

# Load the libraries
library(genefilter)

# Renaming variables
pcut<-p.value.threshold
ecut<-expression.threshold

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=TRUE, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Sanity checks
if(ncol(dat2)==1) {
   stop("You need to have at least two chips to run this analysis!")
}

# Scaling the chip to a mean of 0
dat2<-genescale(dat2)

# Testing
n<-ncol(dat2)
d<-sqrt(n)
s<-rowSds(dat2)/d
m<-rowSums(dat2)/n
T<-(m-0)/s
p.raw<-(1-pt(abs(T), df=n-1))*2
p.adj <- p.adjust(p.raw, method=p.value.adjustment.method)

# Coloring the dot
cols<-rep(1, length(m))
cols[which(p.adj<=pcut & m<(-ecut))] <- 3
cols[which(p.adj<=pcut & m>ecut)] <- 2

# Plotting
pdf(file="volcanoP.pdf", width=w/72, height=h/72)
plot(m, -log10(p.adj), xlim=c(-max(m), max(m)), ylim=c(0,ceiling (max(-log10(p.adj)))), col=cols, main="Volcano plot", pch=19, xlab="Mean expression", ylab="-log10 (p)")
dev.off()

pdf(file="volcanoSE.pdf", width=w/72, height=h/72)
symbols(m, -log10(p.adj), rectangles=cbind(rep(0, length(p.adj)), s), fg=cols, xlim=c(-max(m), max(m)), ylim=c(0,ceiling (max(-log10(p.adj)))), main="Volcano plot", xlab="Mean expression", ylab="-log10 (p)")
dev.off()

# Printing output table
write.table(data.frame(dat, p.adjusted=round(p.adj, digits=6), FC=m), file="one-sample.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

