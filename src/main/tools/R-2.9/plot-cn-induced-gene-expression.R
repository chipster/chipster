# ANALYSIS "aCGH tools (beta testing)"/"Plot copy-number-induced gene expression" (Plot the expression levels of an individual gene for a copy number vs. expression comparison (using the cn-induced-expression.RData file).)
# INPUT GENERIC cn-induced-expression.RData
# OUTPUT cn-induced-gene.png
# PARAMETER gene.id INTEGER (The gene.id to be plotted.)

# plot-cn-induced-gene-expression.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-03-11

library(intCNGEan)

load('cn-induced-expression.RData')

bitmap(file='cn-induced-gene.png', width=600/72, height=600/72)
intCNGEan.plot(gene.id=gene.id, result, tuned)
dev.off()

# EOF