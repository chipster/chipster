# ANALYSIS Statistics/"Sample size estimation" (Estimates sample size on the basis of a set of normalized control 
# samples.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT size.pdf, power.pdf, delta.pdf
# PARAMETER effect.size DECIMAL FROM 0 TO 100 DEFAULT 2 (Fold change)
# PARAMETER alpha DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value) 
# PARAMETER beeta DECIMAL FROM 0 TO 1 DEFAULT 0.8 (Power)
# PARAMETER group.size INTEGER FROM 1 TO 1000 DEFAULT 4 (Sample size)


# JTT 15.12.2006

# Parameter settings (default) for testing purposes
#effect.size<-c(2)
#alpha<-c(0.05)
#beeta<-c(0.80)
#group.size<-c(4)

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

bitmap(file="size.pdf", width=600/72, height=600/72)
ssize.plot(size, xlim=c(0,20), main=paste("Sample size to detect ", effect.size, "-fold change", sep=""))
dev.off()

bitmap(file="power.pdf", width=600/72, height=600/72)
power.plot(po, main=paste("Power to detect ", effect.size, "-fold change", sep=""))
dev.off()

bitmap(file="delta.pdf", width=600/72, height=600/72)
delta.plot(del, xlim=c(0,20), main=paste("Fold change to achieve ", beeta*100, "% power", sep=""))
dev.off()
