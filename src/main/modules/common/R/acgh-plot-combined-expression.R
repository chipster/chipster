# TOOL acgh-plot-combined-expression.R: "Plot copy-number-induced gene expression" (Plot the expression levels of individual genes for a copy number vs. expression comparison. This tool must be run on the output from the tool Test for copy number induced expression changes.)
# INPUT cn-induced-expression.tsv: cn-induced-expression.tsv TYPE GENE_EXPRS 
# OUTPUT cn-induced-expression-plot.pdf: cn-induced-expression-plot.pdf 
# PARAMETER gene.ids: "Gene IDs" TYPE STRING DEFAULT 1 (The gene.ids of the genes to be plotted, separated by commas. Ranges are also supported (e.g. 1,3,7-10\). Ignored, if field genes is not empty.)
# PARAMETER genes: "Gene symbols" TYPE STRING DEFAULT empty (Gene symbols or probe names to be plotted, separated by commas. If this field is filled, the gene.ids parameter will be ignored.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2013-03-25

source(file.path(chipster.common.path, "library-Chipster.R"))
source(file.path(chipster.common.path, "library-intCNGEan.R"))

# read input file
input <- readData("cn-induced-expression.tsv")

tuned <- toIntCNGEanTuned(input)

# parse what to plot
if (genes == "" || genes == "empty") {
  to.plot <- parseSamplesToPlot(gene.ids, input$gene.id)
} else {
  items <- strsplit(genes, ",", genes)[[1]]
  to.plot <- integer()
  for (item in items) {
    to.plot <- c(to.plot, input[gsub(" ", "", item), "gene.id"])
    to.plot <- c(to.plot, input[input$symbol == gsub(" ", "", item), "gene.id"])
  }
  to.plot <- unique(to.plot)
}

# check that we have something to plot
if (length(to.plot) == 0)
  stop("CHIPSTER-NOTE: Nothing to plot.")

# plot
pdf(file="cn-induced-expression-plot.pdf", paper="a4r", width=0, height=0)
for (gene in to.plot)
  intCNGEan.plot(gene.id=gene, tuned)
dev.off()

# EOF
