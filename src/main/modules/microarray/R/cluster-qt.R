# TOOL cluster-qt.R: "Quality Threshold (QT\)" (Quality threshold clustering for genes. Divides the genes of the selected data set in a number of clusters. The number and tightness of the clusters is controlled by the option radius. A suitable radius is dependent on the variability of the data set, and it may vary from one data set to another.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT qt.tsv: qt.tsv 
# OUTPUT qt.pdf: qt.pdf 
# PARAMETER radius.for.similarity: "Radius for similarity" TYPE DECIMAL FROM 0 TO 1000 DEFAULT 20 (A radius for similar genes)
# PARAMETER image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the resampling image)
# PARAMETER image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the resampling image)


# K-means clustering
# JTT 22.6.2006

# Parameter settings (default) for testing purposes
#radius.for.similarity<-20
#image.width<-600
#image.height<-600

# Renaming variables
rad<-radius.for.similarity
w<-image.width
h<-image.height

# Load the libraries
library(flexclust)

# Loads the data file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]
dat3<-dat2

# Calculates the K-means clustering result
# Needs a parameter for groups radius
# rad<-c(20)
qtc<-try(qtclust(dat3, radius=rad))
if(class(qtc) == "try-error") {
	if(length(grep("one cluster", qtc[1]))==1) {
		stop("CHIPSTER-NOTE: All points in one cluster, please choose a smaller radius")
	}
	if(length(grep("Could not find a valid clustering", qtc[1]))==1) {
		stop("CHIPSTER-NOTE: Too many clusters, please choose a higher radius")
	}
}

# Plotting the clustering
max.dat2<-max(dat3)
min.dat2<-min(dat3)
k<-length(unique(qtc@cluster))

counter <- 0
a <- try(log("a"))
while(counter < 5 && class(a) == "try-error") {
	w.new <- w + (counter*100);
	h.new <- h + (counter*100);
	
	pdf(file="qt.pdf", width=w.new/72, height=h.new/72)
	par(mfrow=c(ceiling(sqrt(k)), ceiling(sqrt(k))))

	for(i in 1:k) {
   		a <- try(matplot(t(dat3[qtc@cluster==i,]), type="l", main=paste("cluster:", i), ylab="log expression", col=1, lty=1, ylim=c(min.dat2, max.dat2)))
	}
	dev.off()
	counter <- counter + 1;
}

if(class(a)=="try-error") {
	stop(paste("CHIPSTER-NOTE: Could not draw cluster images even thought width was set to", w.new, "and height was set to", h.new,". Please make image area larger", sep=" "))
}

# Writing a table
# Creates a table with one column giving the cluster membership
dat4<-data.frame(dat3,cluster=qtc@cluster)
dat4<-na.omit(dat4)
write.table(dat4, "qt.tsv", sep="\t", row.names=T, col.names=T, quote=F)
