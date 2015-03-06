# TOOL impute.R: "Impute missing values" (Imputation of missing values. If the maximum specified number of missing values is exceeded, the missing values for that stratum of data are not replaced with imputed values.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT imputed.tsv: imputed.tsv 
# PARAMETER imputation.method: "Imputation method" TYPE [mean: mean, median: median, knn: knn] DEFAULT knn (Imputation method)
# PARAMETER number.of.neighbors: "Number of neighbors" TYPE INTEGER FROM 1 TO 100000 DEFAULT 5 (Number of neighbors to use for knn imputation)
# PARAMETER missing.row.max: "Missing values row max" TYPE PERCENT DEFAULT 20 (Maximum number of missing values on a row)
# PARAMETER missing.column.max: "Missing values column max" TYPE PERCENT DEFAULT 20 (Maximum number of missing values on a column)


# JTT 22.06.2006: Imputation of missing values by mean or median
# IS  12.10.2012: to cope with tables with gene descriptions (that typically contain 's)

# Parameter settings (default) for testing purposes
#imputation.method<-c("mean")
#number.of.neighbors<-c(5)
#missing.row.max<-c(20)
#missing.column.max<-c(20)

# Renaming variables
method<-imputation.method
K<-number.of.neighbors
rmax<-missing.row.max/100
cmax<-missing.column.max/100

# Loads the file
file <- 'normalized.tsv'
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, check.names=FALSE)

# Separates expression values and flags
calls <- dat[,grep('^flag', names(dat))]
dat2 <- dat[,grep('^chip', names(dat))]
# annotations<-dat[,grep("annotations", names(dat))]
annotations <- dat[,grep('^flag|^chip|^average', names(dat), invert=TRUE)]
A <- dat[,grep('^average', names(dat))]

# Inf/-Inf values replaced by NAs
dat2[is.infinite(as.matrix(dat2))]<-NA

# Imputation
if(method=="mean" | method=="median") {
	library(e1071)
	dat.impute<-impute(dat2, method)
	dat.impute<-data.frame(dat.impute, calls)
}
if(method=="knn") {
	library(impute)

	if(max(apply(dat2, 2, function(z) sum(is.na(z))) / nrow(dat2)) > cmax) {
		namax <- max(apply(dat2, 2, function(z) sum(is.na(z))) / nrow(dat2))
		stop(paste("CHIPSTER-NOTE: One of your columns has more NAs than anticipated. Please choose another imputation method or set missing.row.max to ", round(namax, 2)*100, "%", sep=""))
	}

	dat.impute<-impute.knn(as.matrix(dat2), k = K, rowmax = rmax, colmax = cmax, maxp = "p")
	dat.impute<-data.frame(dat.impute$data, calls)
}
if(ncol(A)>0) {
	dat.impute<-data.frame(dat.impute, A)
}

# Add the annotations and flags back into the data table
dat3 <- cbind(annotations,dat.impute)
dat4 <- cbind(dat3,calls)

# Writes a table 
options(scipen=10)
write.table(dat4, file=("imputed.tsv"), sep="\t", row.names=T, col.names=T, quote=F)

# EOF
