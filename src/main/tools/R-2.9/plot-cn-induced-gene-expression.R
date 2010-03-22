# ANALYSIS "aCGH tools (beta testing)"/"Plot copy-number-induced gene expression" (Plot the expression levels of individual genes for a copy number vs. expression comparison (using the cn-induced-expression.RData file).)
# INPUT GENERIC cn-induced-expression.RData
# OUTPUT cn-induced-gene.png
# PARAMETER genes STRING DEFAULT 1 (The gene.ids of the genes to be plotted, separated by commas. Ranges are also supported (e.g. 1,3,7-10).)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# plot-cn-induced-gene-expression.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-03-19

library(intCNGEan)

# read input file
load('cn-induced-expression.RData')

# parse the input string
genes <- gsub('[^0-9,-]', ',', genes)
items <- strsplit(genes, ',')[[1]]
to.plot <- integer()
for (item in items) {
	item <- item[item!='']
	if (length(item)==0) next
	range <- strsplit(item, '-')[[1]]
	range <- range[range!='']
	if (length(range)==0) next
	to.plot <- c(to.plot, seq(range[1], range[length(range)]))
}
to.plot <- unique(to.plot)

# remove genes that are out of bounds
to.plot <- to.plot[to.plot<=max(result$gene.id)]

# check that we have something to plot
if (length(to.plot)==0)
	stop('Nothing to plot.')

# plot
bitmap(file='cn-induced-gene.png', width=image.width/72, height=image.height/72)
if (length(to.plot)==1) {
	intCNGEan.plot(gene.id=to.plot, result, tuned)
} else {
	sq <- sqrt(length(to.plot))
	rows <- ceiling(sq)
	cols <- ceiling(length(to.plot)/rows)
	par(mfrow=c(rows,cols))
	for (gene in to.plot)
		intCNGEan.plot(gene.id=gene, result, tuned)
}
dev.off()

# EOF