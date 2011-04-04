# ANALYSIS Preprocessing/"Filter by CV" (Filter genes by their coefficient of variation.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT cv-filter.tsv
# PARAMETER percentage.to.filter.out DECIMAL FROM 0 TO 1 DEFAULT 0.5 (Percentage to filter out)

# JTT, 14.1.2008
#
# MG, 27.12.2010, modified to handle data with zero mean

# Renaming variables
percentage<-percentage.to.filter.out

# Loading the libraries
library(genefilter)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]
calls<-data.frame(calls)
dat2<-data.frame(dat2)

# Sanity checks
if(ncol(dat2)==1) {
   stop("You need to have at least two chips to apply filtering by genes!")
}

# Filter
colno<-ncol(dat2)
g.mean <- rowSums(dat2)/colno
g.sd <- rowSds(dat2)
g.cv <- g.sd / g.mean
	# if an infinite result is obtained because of a zero mean value
	# force the CV value to become zero
	for (count in 1:length(g.cv)) {
		if (is.na (g.cv[count])) {
			g.cv[count] <- 0
		}
	}
sel<-(g.cv>quantile(g.cv,percentage,na.rm=TRUE))
set<-dat[sel, ]

# Saving the results
write.table(data.frame(set), file=("cv-filter.tsv"), sep="\t", row.names=T, col.names=T, quote=F)

