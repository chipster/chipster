# TOOL plot-annoheatmap.R: "Annotated heatmap" (Draws an annotated heatmap. Clustering is done using the function hclust.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT heatmap.pdf: heatmap.pdf
# PARAMETER column: "Annotation column" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column to be used for the annotation)
# PARAMETER number.of.groups: "Number of groups" TYPE INTEGER FROM 0 TO 20 DEFAULT 2 (How many groups to color to the tree)
# PARAMETER OPTIONAL coloring.scheme: "Coloring scheme" TYPE [Green-Red: Green-Red, Green-Black-Red: Green-Black-Red, Blue-White-Red: Blue-White-Red, Black-White: Black-White, None: None] DEFAULT Green-Red (Coloring scheme for the heatmap. Set to None to remove the heatmap entirely from the plot.)
# PARAMETER OPTIONAL cluster.samples.only: "Cluster samples only" TYPE [yes: yes, no: no] DEFAULT yes (Disables clustering on the genes. This option is convenient if you want to retain a predefined gene order or make a sample clustering heatmap with more than 10000 genes.)
# PARAMETER OPTIONAL hm.scale: "Scale data" TYPE [none: none, row: row, column: column] DEFAULT row (Indicates if the values should be centered and scaled in either the row direction or the column direction, or none. Affects only data visualistion, not the actual clustering.)
# PARAMETER OPTIONAL distance: "Distance" TYPE [euclidean: Euclidean,  manhattan: Manhattan, binary: Binary, pearson: Pearson, spearman: Spearman, kendall: Kendall] DEFAULT pearson (The correlation distance measure to be used for clustering.)
# PARAMETER OPTIONAL clu.method: "Clustering method" TYPE [ward: Ward, single: single, complete: complete, average: average, median: median] DEFAULT average (The agglomeration to be used for clustering.)
# PARAMETER OPTIONAL image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the image.)
# PARAMETER OPTIONAL image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the image.)

# MK 17.10.2013 Annotated Heatmap
library("Heatplus")

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Sanity checks
if (nrow(dat2) > 20000 & cluster.samples.only=="no") {
	stop("CHIPSTER-NOTE: Hierarchical clustering can be run on a maximum of 20000 genes");
}
if (ncol(dat2) > 20000) {
	stop("CHIPSTER-NOTE: Hierarchical clustering can be run on a maximum of 20000 samples");
}

# Read phenodata and generate annotation features
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
colnames(dat2)<-gsub(" ", "", phenodata$description)
if(column != "EMPTY") {
	groups<-phenodata[,pmatch(column,colnames(phenodata))]
	mean.dat <- round(apply(dat2, 2, mean),2)
	#feat <- setNames(data.frame(data1=groups, data2=mean.dat), c(column, "mean"))
	feat <- setNames(data.frame(data1=groups), column)
	feat <- convAnnData(feat, nval.fac=ncol(dat2)-1)
	feat <- list(asIs=TRUE, Col=list(data=feat))

	#feat <- list(Col = list(data=feat, fun=picketPlot))
} else {
	feat <- NULL
}

# Clutering parameters and clustering
if(distance=="euclidean" || distance=="manhattan" || distance=="binary") {
	corrdist = function(x) as.dist(dist(x, method=distance))
} else {
	corrdist = function(x) as.dist(1-cor(t(x), method=distance))
}

hclust.met = function(x) hclust(x, method=clu.method)

if (cluster.samples.only=="yes" || coloring.scheme == "None") {
	dend.met = list(Col=list(status="yes", clustfun=hclust.met, distfun=corrdist), Row=list(status="no", clustfun=hclust.met, distfun=corrdist))
} else {
	dend.met = list(Col=list(status="yes", clustfun=hclust.met, distfun=corrdist), Row=list(status="yes", clustfun=hclust.met, distfun=corrdist))
}	

# Tree cutting
if(number.of.groups > 1) {
	ngroups <- number.of.groups - 1
	hc.temp <- hclust(corrdist(t(dat2)), method=clu.method)
	cuth.temp <- hc.temp$height[length(hc.temp$height)-ngroups]
	cuth.info <- list(Col=list(cuth=cuth.temp))
} else {
	cuth.info <- NULL
}

# Data scaling. Used only for visualization purposes
if (hm.scale == "row") {
	dat3 <- dat2
	dat3 <- sweep(dat3, 1, rowMeans(dat3, na.rm = TRUE))
	sd <- apply(dat3, 1, sd, na.rm = TRUE)
	dat3 <-sweep(dat3, 1, sd, "/")
	dat3[dat3 > 3] <- 3;
	dat3[dat3 < -3] <- -3;
	plot.range <- c(-3,3);
} else if (hm.scale == "column") {
	dat3 <- dat2
	dat3 <- sweep(dat3, 2, colMeans(dat3, na.rm = TRUE))
	sd <- apply(dat3, 2, sd, na.rm = TRUE)
	dat3 = sweep(dat3, 2, sd, "/")
	dat3[dat3 > 3] <- 3;
	dat3[dat3 < -3] <- -3;
	plot.range <- c(-3,3);
} else {
	dat3 <- dat2
	plot.range <- c(min(dat2),max(dat2));
}

# Color scheme
num.breaks <- 512
break.size <- length(niceBreaks(plot.range, num.breaks))-1

if(coloring.scheme == "Green-Red" || coloring.scheme == "None") {
	colmap <- colorRampPalette(c("Green", "Red"), space = "rgb")(break.size)
} else if (coloring.scheme == "Green-Black-Red") {
	colmap<-colorRampPalette(c("Green", "Black", "Red"))(break.size)
} else if (coloring.scheme == "Blue-White-Red") {
	colmap<-colorRampPalette(c("Blue", "White", "Red"))(break.size)
} else if (coloring.scheme == "Black-White") {
	colmap<-colorRampPalette(c("Black", "LightGrey"))(break.size)
}

# Set column spacing
longest.name <- colnames(dat2)[which.max(unlist(lapply(colnames(dat2), nchar)))];
temp.width <- floor(strwidth(longest.name, unit="in", cex=(0.2 + 1/log10(ncol(dat2)))) / par("csi") + 1)

# Pseudo-plot to find par("csi") value
pdf(file="heatmap.pdf", width=image.width/72, height=image.height/72)
ann1 = annHeatmap2(as.matrix(dat2), col = colmap, breaks = niceBreaks(plot.range, num.breaks), scale="none", dendrogram=dend.met, annotation=feat, cluster=cuth.info, labels=list(Col=list(nrow=temp.width + 1)), legend=T)
try_plot <- try(plot(ann1), silent=T)
if(class(try_plot) == "try-error") {
	stop("CHIPSTER-NOTE: Your plot area is too small. Please consider increasing the width and height of the plot")
}
longest.name <- colnames(dat2)[which.max(unlist(lapply(colnames(dat2), nchar)))];
temp.width <- floor(strwidth(longest.name, unit="in", cex=(0.2 + 1/log10(ncol(dat2)))) / par("csi") + 1)
dev.off()

# Create plot
if(coloring.scheme != "None") {
	pdf(file="heatmap.pdf", width=image.width/72, height=image.height/72)
	ann1 = annHeatmap2(as.matrix(dat2), col = colmap, breaks = niceBreaks(plot.range, num.breaks), scale="none", dendrogram=dend.met, ann=feat, cluster=cuth.info, labels=list(Col=list(nrow=temp.width + 1)), legend=T)
	#done like this, as no idea how to otherwise visualise values beyond plot.range
	ann1$data$x2 <- dat3[match(rownames(ann1$data$x2), rownames(dat3)), match(colnames(ann1$data$x2), colnames(dat3))]
	plot(ann1)
	dev.off()
} else {
	#values above and below plot.range are white
	pdf(file="heatmap.pdf", width=image.width/72, height=image.height/72)
	#ann1 = annHeatmap2(as.matrix(dat2), col = colmap, breaks = niceBreaks(plot.range, num.breaks), scale="none", dendrogram=dend.met, annotation=feat, cluster=cuth.info, labels=list(Col=list(nrow=temp.width + 1), Row=list(labels=rep(" ", nrow(dat2)))), legend=F)
	ann1 = annHeatmap2(as.matrix(dat2), col = colmap, breaks = niceBreaks(plot.range, num.breaks), scale="none", dendrogram=dend.met, ann=feat, cluster=cuth.info, labels=list(Col=list(nrow=temp.width + 1), Row=list(labels=rep(" ", nrow(dat2)))), legend=F)

	ann1$data$x2[ann1$data$x2 < max(plot.range)] <- max(plot.range) + 1
	plot(ann1, widths=c(2,5,1), heights=c(2,0.75,1))
	dev.off()
}
