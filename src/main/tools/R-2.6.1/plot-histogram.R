# ANALYSIS Visualisation/"Histogram" (Creates a histogram for every chip using normalized data.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT histogram.png
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Histogram
# JTT 4.10.2007

# Renaming variables
w<-image.width
h<-image.height

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Draws the histograms
s<-ceiling(sqrt(ncol(dat2)))
titles<-gsub(" ", "", phenodata$description)

if(length(titles)==ncol(dat2)) {
   bitmap(file="histogram.png", width=w/72, height=h/72)
   par(mfrow=c(s,s))
   for(i in 1:ncol(dat2)) {
      par(mar=c(2,2,3,0))
      hist(dat2[,i], main=titles[i], cex.main=0.80, br=100, col=1, border=1)
   }
   dev.off()
} else {
   bitmap(file="histogram.png", width=w/72, height=h/72)
   par(mfrow=c(s,s))
   for(i in 1:ncol(dat2)) {
      par(mar=c(2,2,3,0))
      hist(dat2[,i], cex.main=0.80, br=100, col=1, border=1)
   }
   dev.off()
}
