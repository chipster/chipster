# ANALYSIS Preprocessing/"Filter by standard deviation" (Filtering the genes according to their standard deviation. 
# Specify the cut-off in percents. Deviation is expressed as a percentage of data to filter out. These can be made to
# approximately correspond to SDs assuming normality of the data: 1SD=67%, 2SDs=95%, 3SDs=99.7%.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT sd-filter.tsv
# PARAMETER percentage.to.filter.out DECIMAL FROM 0 TO 1 DEFAULT 0.997 (Percent of values to filter out)

# JTT, 13.7.2005
# Heavily modified on 7.6.2006
# Even heavier modifications on 20.9.2007

# Renaming variables
percentage<-percentage.to.filter.out

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

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
write.table(data.frame(set), file=("sd-filter.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
