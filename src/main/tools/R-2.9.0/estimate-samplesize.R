# ANALYSIS Statistics/"Sample size estimation" (Estimates sample size on the basis of a set of normalized control 
# samples.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT size.png, power.png, delta.png
# PARAMETER effect.size DECIMAL FROM 0 TO 100 DEFAULT 2 (Fold change)
# PARAMETER alpha DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value) 
# PARAMETER beeta DECIMAL FROM 0 TO 1 DEFAULT 0.8 (Power)
# PARAMETER group.size INTEGER FROM 1 TO 1000 DEFAULT 4 (Sample size)


# JTT 15.12.2006

# Loads the libraries
library(genefilter)
library(ssize)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Calculate standard deviation row-wise
sds<-rowSds(dat2)

# Calculate the sample size, power and delta
size<-ssize(sd=sds, delta=log2(effect.size), sig.level=alpha, power=beeta)
po<-pow(sd=sds, n=group.size, delta=log2(effect.size), sig.level=alpha)
del<-delta(sd=sds, n=group.size, power=beeta, sig.level=alpha)

# png("size.png", width=600, height=600)
bitmap(file="size.png", width=600/72, height=600/72)
   ssize.plot(size, xlim=c(0,20), main=paste("Sample size to detect ", effect.size, "-fold change", sep=""))
dev.off()

# png("power.png", width=600, height=600)
bitmap(file="power.png", width=600/72, height=600/72)
   power.plot(po, main=paste("Power to detect ", effect.size, "-fold change", sep=""))
dev.off()

# png("delta.png", width=600, height=600)
bitmap(file="delta.png", width=600/72, height=600/72)
   delta.plot(del, xlim=c(0,20), main=paste("Fold change to achieve ", beeta*100, "% power", sep=""))
dev.off()
