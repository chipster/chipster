# ANALYSIS "aCGH tools (beta testing)"/"Plot combined profiles of copy number and expression" (Plot the copy number vs. expression profile of an individual array (using the cn-induced-expression.RData file).)
# INPUT GENERIC cn-induced-expression.RData
# OUTPUT cn-induced-profile.png
# PARAMETER sample STRING DEFAULT 1 (The number of the sample to be plotted.)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# plot-cn-induced-expression-profile.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-03-20

# The plotting command from the intCNGEan package uses par(mfrow) overriding ours.
# Therefore only one profile can be plotted.
# However, the code already contains the necessary functionality to deal with multiple samples,
# in case I find a way how to do it in the future.
# If that happens, the parameter definition will be altered, and this line removed:
samples <- sample

library(intCNGEan)

# read input file
load('cn-induced-expression.RData')

# parse the input string
samples <- gsub('[^0-9,-]', ',', samples)
items <- strsplit(samples, ',')[[1]]
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

# remove samples that are out of bounds
to.plot <- to.plot[to.plot<=length(sampleNames(matched$CNdata.matched))]

# check that we have something to plot
if (length(to.plot)==0)
	stop('Nothing to plot.')

# plot
bitmap(file='cn-induced-profile.png', width=image.width/72, height=image.height/72)
if (length(to.plot)==1) {
	intCNGEan.profilesPlot(matched$CNdata.matched, matched$GEdata.matched, sampleNo=to.plot)
} else {
	sq <- sqrt(length(to.plot))
	rows <- ceiling(sq)
	cols <- ceiling(length(to.plot)/rows)
	par(mfrow=c(rows,cols))
	for (sample in to.plot)
		intCNGEan.profilesPlot(matched$CNdata.matched, matched$GEdata.matched, sampleNo=sample)
}
dev.off()

# EOF