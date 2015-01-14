# TOOL filter-sd.R: "Filter by standard deviation" (Filters genes according to their standard deviation. Specify what percentage of the genes showing the lowest standard deviation should be filtered out. These can be made to approximately correspond to SDs assuming normality of the data: 1SD=67%, 2SDs=95%, 3SDs=99.7%. The default is 50%.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT sd-filter.tsv: sd-filter.tsv 
# PARAMETER percentage.to.filter.out: "Percentage to filter out" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.50 (Percent of genes showing lowest standard deviation to be filtered out.)

# JTT, 13.7.2005
# Heavily modified on 7.6.2006
# Even heavier modifications on 20.9.2007
# IS, 16.10.2012, modified to cope with tables with gene descriptions (that typically contain 's)
# EK, 5.11.2014, changed the default value and clarified the text.

# Parameter settings (default) for testing purposes
#percentage.to.filter.out<-c(0.997)

# Renaming variables
percentage<-percentage.to.filter.out

# Loads the normalized data
file <- c('normalized.tsv')
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, check.names=FALSE)

# Sanity checks
if(ncol(dat)==1) {
   stop("You need to have at least two chips to apply filtering by genes!")
}

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]
calls<-data.frame(calls)
dat2<-data.frame(dat2)

# Filter
sds<-apply(as.data.frame(dat2), MARGIN=1, FUN="sd")
sel<-(sds>quantile(sds,percentage))
set<-dat[sel, ]

# Saving the results
options(scipen=10)
write.table(data.frame(set), file=("sd-filter.tsv"), sep="\t", row.names=T, col.names=T, quote=F)

# EOF
