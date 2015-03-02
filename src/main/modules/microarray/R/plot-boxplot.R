# TOOL plot-boxplot.R: Boxplot (Creates a boxplot of normalized data. One box per chip is plotted.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT boxplot.png: boxplot.png 
# PARAMETER column: Column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column according to which boxplots are coloured)
# PARAMETER OPTIONAL image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER OPTIONAL image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# JTT 02.10.2007: Boxplot
# MG 15.09.2011: Updated colors and legend
# MK 10.10.2013: User can choose which phenodata column is used for coloring bars

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
if(length(grep(column, colnames(phenodata))) == 0) {
	stop("CHIPSTER-NOTE: You need to select phenodata column to run this analysis")
}
phenodata_col <- phenodata[, which(column==colnames(phenodata))]

# Setup sample colors according to group
if (length(levels(as.factor(phenodata_col))) > 0) {
	sample_colors <- numeric(length(phenodata_col))
	group_levels <- levels(as.factor(phenodata_col))
	group_identity <- as.character(phenodata_col)
	for (count_levels in 1:length(group_levels)) {
		for (count_samples in 1:length(phenodata_col)) {
			if(group_identity[count_samples]==group_levels[count_levels]) sample_colors[count_samples] <- 1+count_levels
		}
	}
level_colors <- levels(as.factor(sample_colors))
}
if (length(levels(as.factor(phenodata_col))) == 0) {
	sample_colors <- rep(2,length(phenodata_col))
}

# Plotting
if(nrow(phenodata)==ncol(dat2)) {
	boxplot_names = phenodata$description
} else {
	boxplot_names = colnames(dat2)
}

bitmap(file="boxplot.png", width=w/72, height=h/72)
par(mar=c(12,5,5,5))
boxplot(as.data.frame(dat2), las=2, names=boxplot_names, col=sample_colors)
if (length(levels(as.factor(phenodata_col))) > 0) {
	legend (x="topleft", legend=group_levels, col=level_colors, cex=0.5, pch=19)
}
dev.off()

