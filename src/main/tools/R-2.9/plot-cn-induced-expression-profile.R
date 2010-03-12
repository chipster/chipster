# ANALYSIS "aCGH tools (beta testing)"/"Plot profile of copy-number-induced expression" (Plot the copy number vs. expression profile of an individual array (using the cn-induced-expression.RData file).)
# INPUT GENERIC cn-induced-expression.RData
# OUTPUT cn-induced-profile.png
# PARAMETER sample INTEGER (The number of the sample to be plotted.)

# plot-cn-induced-expression-profile.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-03-11

library(intCNGEan)

load('cn-induced-expression.RData')

bitmap(file='cn-induced-profile.png', width=600/72, height=600/72)
intCNGEan.profilesPlot(matched$CNdata.matched, matched$GEdata.matched, sampleNo=sample)
dev.off()

# EOF