# TOOL cluster-som.R: "Self-organizing map (SOM\)" (Self-organizing map clustering of genes. Note that you need to have at least two chips, and columns times rows number of genes in order to able to execute the analysis.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT som.tsv: som.tsv 
# OUTPUT som.pdf: som.pdf 
# PARAMETER number.of.columns: number.of.columns TYPE INTEGER FROM 1 TO 100 DEFAULT 4 (Number of SOM grid columns)
# PARAMETER number.of.rows: number.of.rows TYPE INTEGER FROM 1 TO 100 DEFAULT 3 (Number of SOM grid rows)
# PARAMETER coloring.scheme: coloring.scheme TYPE [Red-Green: Red-Green, Blue-Yellow: Blue-Yellow, Black-White: Black-White] DEFAULT Blue-Yellow (Coloring scheme for the SOM map)
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the resampling image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the resampling image)


# Self-organizing map
# JTT 22.6.2006

# Parameter settings (default) for testing purposes
#number.of.columns<-4
#number.of.rows<-3
#coloring.scheme<-c("Blue-Yellow")
#image.width<-600
#image.height<-600

# Renaming variables
x<-number.of.columns
y<-number.of.rows
colpar<-coloring.scheme
w<-image.width
h<-image.height

# Loads the data file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Loads the libraries
library(kohonen)
library(RColorBrewer)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Sanity checks
if(x*y>=nrow(dat2)) {
   stop("The number of nodes is larger than the number of genes! Please use a smaller grid.")
}
if(ncol(dat2)<2) {
   stop("You have have at least two chips to run this analysis!")
}

# Calculating the SOM-result
iter<-c(10000)
top<-c("hexagonal")
kohmap<-som(as.matrix(dat2), grid=somgrid(xdim=x, ydim=y, topo=top), rlen=iter)
dists <- unit.distances(kohmap$grid, toroidal=F)

# Makes colors for plotting purposes
dif.dists<-unique(as.vector(names(table(dists[1,]))))
if(colpar=="Red-Green") {
   heatcol<-colorRampPalette(c("Red", "Green"))(length(dif.dists))
}
if(colpar=="Blue-Yellow") {
   heatcol<-colorRampPalette(c("Blue", "Yellow"))(length(dif.dists))
}
if(colpar=="Black-White") {
   heatcol<-colorRampPalette(c("Black", "LightGrey"))(length(dif.dists))
}
heatcol2<-dists[1,]
for(i in 1:length(dif.dists)) {
   heatcol2[heatcol2==dif.dists[i]]<-heatcol[i]
}

# Plotting the SOM map
pdf(file="som.pdf", width=w/72, height=h/72)
par(mfrow=c(2,2), mar=c(2,2,0,0))
plot.kohonen(kohmap, type="property", property=dists[1,], palette.name=heat.colors, main="Property")
plot.kohonen(kohmap, type="codes", property=dists[1,], palette.name=heat.colors, main="Codes")
plot.kohonen(kohmap, type="counts", property=dists[1,], palette.name=heat.colors, main="Counts")
plot.kohonen(kohmap, type="mapping", property=dists[1,], palette.name=heat.colors, pch=19, main="Mapping")
dev.off()

# Makes a dataframe of data, SOM result (cell membership and distance to first cell), and colors
# Dists, cols, and griddim are not of data length, ends are padded with white spaces (" ")
dists<-c(dists[1,], rep(" ", (dim(dat)[[1]]-length(dists[1,]))))
cols<-c(heatcol2, rep(" ", (dim(dat)[[1]]-length(heatcol2))))
griddim<-c(x, y, rep(" ", (dim(dat)[[1]]-2)))
write.table(data.frame(dat, cluster=kohmap$unit.classif, distance2first=dists, colours=cols, griddim=griddim), "som.tsv", sep="\t", row.names=T, col.names=T, quote=F)

