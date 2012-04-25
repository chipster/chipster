# TOOL plot-boxplot.R: Boxplot (Creates a boxplot of normalized data. One box per chip is plotted.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT boxplot.png: boxplot.png 
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# Boxplot
# JTT 2.10.2007
# MG, 15.9.2011
# Updated colors and legend

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

# Setup sample colors according to group
if (length(levels(as.factor(phenodata$group))) > 0) {
	sample_colors <- numeric(length(phenodata$group))
	group_levels <- levels(as.factor(phenodata$group))
	group_identity <- as.character(phenodata$group)
	for (count_levels in 1:length(group_levels)) {
		for (count_samples in 1:length(phenodata$group)) {
			if(group_identity[count_samples]==group_levels[count_levels]) sample_colors[count_samples] <- 1+count_levels
		}
	}
level_colors <- levels(as.factor(sample_colors))
}
if (length(levels(as.factor(phenodata$group))) == 0) {
	sample_colors <- rep(2,length(phenodata$group))
}

# Plotting
if(nrow(phenodata)==ncol(dat2)) {
	bitmap(file="boxplot.png", width=w/72, height=h/72)
	par(mar=c(12,5,5,5))
	boxplot(as.data.frame(dat2), las=2, names=phenodata$description, col=sample_colors)
	if (length(levels(as.factor(phenodata$group))) > 0) {
		legend (x="topleft", legend=group_levels, col=level_colors, cex=0.5, pch=19)
	}
	dev.off()
} else {
	bitmap(file="boxplot.png", width=w/72, height=h/72)
	par(mar=c(12,5,5,5))
	boxplot(as.data.frame(dat2), las=2, col=sample_colors)
	if (length(levels(as.factor(phenodata$group))) > 0) {
		legend (x="topleft", legend=group_levels, col=level_colors, cex=0.5, pch=19)
	}
	dev.off()
}

