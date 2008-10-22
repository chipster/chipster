# ANALYSIS Pathways/"*Bayesian network" (Infers a bayesian network from timeseries. Specify a p-value cut-off for 
# significance, and width and height of the network image. You need to have a time column in your phenodata, 
# otherwise the analysis will fail miserably. Note that this works only, if you have exactly one replicate per 
# each time point.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT network.png
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)
# PARAMETER coloring [solid, expression] DEFAULT solid (Color for nodes in the network)

# Inference of networks
# JTT 15.9.2006
# Corrected 27.9.2006

# Load the libraries
library(e1071)
library(GeneTS)

# Renaming variables
p.cut<-p.value.threshold
w<-image.width
h<-image.height

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Reads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# How many replicates there are per time point?
# This loop needs more work as it currently doesn't handle replicate time points
if(length(phenodata$time)>length(levels(as.factor(phenodata$time)))) {
   repl<-c()
   for(i in 1:length(levels(as.factor(phenodata$time)))) {
      repl<-c(repl, sum(as.numeric(grep(phenodata$time[1], phenodata$time, value=T)==phenodata$time[1])))
   }
} else {
   repl<-rep(1, length(phenodata$time))
}

# Making a longitudinal object
# Corrected to work
dat3<-as.longitudinal(t(dat2), repeats=repl, time=as.numeric(levels(as.factor(phenodata$time))))

# Replacing missing values
dat4<-t(na.omit(t(impute(dat3))))

# Infers and draws the network
# Estimate partial correlation
pcor<-ggm.estimate.pcor(dat4)
# Assign p-values, Q-values and post. probs for edges
pcor.test<-ggm.test.edges(pcor)
# Which genes are significant
significant2.idx <- pcor.test$qval <= p.cut
# Gets gene names
node.labels<-colnames(dat4)
# Generates a graph object using the significant edges
gr <- ggm.make.graph( pcor.test[significant2.idx,], node.labels, drop.singles=TRUE)
# Coloring

# Not yet implemented
 
# Plot the network
png(width=w, height=h, file="network.png")
# bitmap(file="network.png", width=w/72, height=h/72)
ggm.plot.graph(gr,  main="Network", show.edge.labels=TRUE, layoutType=c("neato")) # neato or dot
dev.off()
