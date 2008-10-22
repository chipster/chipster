# ANALYSIS Visualisation/"Volcano plot" (Tests whether the genes are differentially expressed using one-sample t-test,
# and plots the result in a form of a Volcano plot.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT volcanoP.png, volcanoSE.png
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER expression.threshold DECIMAL FROM -10 TO 10 DEFAULT 1 (Expression threshold for plotting) 
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Volcano plot
# JTT 19.9.2007

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
dat<-read.table(file, header=T, sep="\t", row.names=1)

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

# Coloring the dot
cols<-rep(1, length(m))
cols[which(p.raw<=pcut & m<(-ecut))] <- 3
cols[which(p.raw<=pcut & m>ecut)] <- 2

# Plotting
bitmap(file="volcanoP.png", width=w/72, height=h/72)
plot(m, -log(p.raw), xlim=c(-max(m), max(m)), col=cols, main="Volcano plot", pch=19, xlab="Mean expression", ylab="-log(p)")
dev.off()

bitmap(file="volcanoSE.png", width=w/72, height=h/72)
symbols(m, -log(p.raw), rectangles=cbind(rep(0, length(p.raw)), s), ylim=c(0, 40), fg=cols, xlim=c(-max(m), max(m)), main="Volcano plot", xlab="Mean Expression", ylab="-log(p)")
dev.off()


