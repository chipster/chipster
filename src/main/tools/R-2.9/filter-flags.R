# ANALYSIS Preprocessing/"Filter by flags" (Gene filtering using flags. If more than the specified number of
# chips have the specified flag for a particular gene, it passes the filter, and is saved in a new data set.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT call-filtered.tsv
# PARAMETER filter.by.flag [P, M, A, A.or.M, M.or.P] DEFAULT P (Which flag to filter by)
# PARAMETER number.of.chips INTEGER FROM 1 TO 10000 DEFAULT 2 (How many samples the should have the specified call value)


# Filtering using flags (calls)
# JTT 29.6.2006
#
# MG 14.9.2010
# Modified to allow filtering options for A or M and M or P

# Renaming variables
f1 <- filter.by.flag
f2 <- f1
if (filter.by.flag=="A.or.M") {
	f1 <- "A"
	f2 <- "M"
}
if (filter.by.flag=="M.or.P") {
	f1 <- "M"
	f2 <- "P"
}
p<-number.of.chips

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Sanity checks
if(ncol(as.data.frame(calls))==0) {
	stop("You do not have any flags!")
}

# Converting the data to data frames
dat2<-data.frame(dat2)
calls<-data.frame(calls)
len<-length(calls)

# Apply filtering
for(i in 1:len) {
	calls[,i]<-ifelse((calls[,i]==f1 | calls[,i]==f2) , 1, 0);
}
s<-apply(calls, MARGIN=1, FUN="sum")
dat2<-dat[which(s>=p),]

# Saving the results into a text file
write.table(dat2, "call-filtered.tsv", sep="\t", row.names=T, col.names=T, quote=F)
