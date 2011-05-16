# TOOL plot-venn-diagram.R: "Venn diagram" (Draws a Venn diagram for the three selected gene lists. Select the gene lists from the Workflow view of Datasets before running this tool.)
# INPUT normalized1.tsv: normalized1.tsv TYPE GENE_EXPRS 
# INPUT normalized2.tsv: normalized2.tsv TYPE GENE_EXPRS 
# INPUT normalized3.tsv: normalized3.tsv TYPE GENE_EXPRS 
# OUTPUT venn.pdf: venn.pdf 
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Parameter settings (default) for testing purposes
#image.width<-600
#image.height<-600

# Renaming variables
w<-image.width
h<-image.height

# Loads the libraries
library(limma)

# Loads the data
files<-dir()
set1<-read.table(files[1], header=T, sep="\t")
set2<-read.table(files[2], header=T, sep="\t")
set3<-read.table(files[3], header=T, sep="\t")

# Reads just the row names
set1<-unique(rownames(set1))
set2<-unique(rownames(set2))
set3<-unique(rownames(set3))

# Creates the universe
universe <- sort(unique(c(set1, set2, set3)))

# Creates a matrix of terms in every set
Counts <- matrix(0, nrow=length(universe), ncol=3)
colnames(Counts) <- c(files)
for (i in 1:length(universe)) {
   Counts[i,1] <- universe[i] %in% set1
   Counts[i,2] <- universe[i] %in% set2
   Counts[i,3] <- universe[i] %in% set3
}

# The part plotting the universe has been commented off.
# Otherwise the code is taken directly from limma
vennDiagram2<-function (object, include = "both", names, mar = rep(1, 4), cex = 1.5, 
    lwd = 1, circle.col, counts.col, show.include, ...) {
    if (!is(object, "VennCounts")) {
        if (length(include) > 2) 
            stop("Cannot plot Venn diagram for more than 2 sets of counts")
        if (length(include) == 2) 
            object.2 <- vennCounts(object, include = include[2])
        object <- vennCounts(object, include = include[1])
    }
    else if (length(include == 2)) 
        include <- include[1]
    nsets <- ncol(object) - 1
    if (nsets > 3) 
        stop("Can't plot Venn diagram for more than 3 sets")
    if (missing(names)) 
        names <- colnames(object)[1:nsets]
    counts <- object[, "Counts"]
    if (length(include) == 2) 
        counts.2 <- object.2[, "Counts"]
    if (missing(circle.col)) 
        circle.col <- par("col")
    if (length(circle.col) < nsets) 
        circle.col <- rep(circle.col, length.out = nsets)
    if (missing(counts.col)) 
        counts.col <- par("col")
    if (length(counts.col) < length(include)) 
        counts.col <- rep(counts.col, length.out = length(include))
    if (missing(show.include)) 
        show.include <- as.logical(length(include) - 1)
    theta <- 2 * pi * (1:360)/360
    xcentres <- list(0, c(-1, 1), c(-1, 1, 0))[[nsets]]
    ycentres <- list(0, c(0, 0), c(1/sqrt(3), 1/sqrt(3), -2/sqrt(3)))[[nsets]]
    r <- c(1.5, 1.5, 1.5)[nsets]
    xtext <- list(-1.2, c(-1.2, 1.2), c(-1.2, 1.2, 0))[[nsets]]
    ytext <- list(1.8, c(1.8, 1.8), c(2.4, 2.4, -3))[[nsets]]
    old.par <- par(mar = mar)
    on.exit(par(old.par))
    plot(x = 0, y = 0, type = "n", xlim = c(-4, 4), ylim = c(-4, 
        4), xlab = "", ylab = "", axes = FALSE, ...)
    for (circle in 1:nsets) {
        lines(xcentres[circle] + r * cos(theta), ycentres[circle] + 
            r * sin(theta), lwd = lwd, col = circle.col[circle])
        text(xtext[circle], ytext[circle], names[circle], cex = cex)
    }
    switch(nsets, {
        rect(-3, -2.5, 3, 2.5)
        printing <- function(counts, cex, adj, col, leg) {
            text(2.3, -2.1, counts[1], cex = cex, col = col, 
                adj = adj)
            text(0, 0, counts[2], cex = cex, col = col, adj = adj)
            if (show.include) 
                text(-2.3, -2.1, leg, cex = cex, col = col, adj = adj)
        }
    }, {
        rect(-3, -2.5, 3, 2.5)
        printing <- function(counts, cex, adj, col, leg) {
            text(2.3, -2.1, counts[1], cex = cex, col = col, 
                adj = adj)
            text(1.5, 0.1, counts[2], cex = cex, col = col, adj = adj)
            text(-1.5, 0.1, counts[3], cex = cex, col = col, 
                adj = adj)
            text(0, 0.1, counts[4], cex = cex, col = col, adj = adj)
            if (show.include) 
                text(-2.3, -2.1, leg, cex = cex, col = col, adj = adj)
        }
    }, {
        rect(-3, -3.5, 3, 3.3)
        printing <- function(counts, cex, adj, col, leg) {
            #text(2.5, -3, counts[1], cex = cex, col = col, adj = adj)
            text(0, -1.7, counts[2], cex = cex, col = col, adj = adj)
            text(1.5, 1, counts[3], cex = cex, col = col, adj = adj)
            text(0.75, -0.35, counts[4], cex = cex, col = col, 
                adj = adj)
            text(-1.5, 1, counts[5], cex = cex, col = col, adj = adj)
            text(-0.75, -0.35, counts[6], cex = cex, col = col, 
                adj = adj)
            text(0, 0.9, counts[7], cex = cex, col = col, adj = adj)
            text(0, 0, counts[8], cex = cex, col = col, adj = adj)
            if (show.include) 
                text(-2.5, -3, leg, cex = cex, col = col, adj = adj)
        }
    })
    adj <- c(0.5, 0.5)
    if (length(include) == 2) 
        adj <- c(0.5, 0)
    printing(counts, cex, adj, counts.col[1], include[1])
    if (length(include) == 2) 
        printing(counts.2, cex, c(0.5, 1), counts.col[2], include[2])
    invisible()
}

# Plotting the diagram
cols<-c("Red", "Green", "Blue")
pdf(file="venn.pdf", width=w/72, height=h/72)
vennDiagram2(vennCounts(Counts), circle.col=cols)
dev.off()
