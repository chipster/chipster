# ANALYSIS Preprocessing/"Impute missing values" (Imputation of missing values. If the maximum specified number of 
# missing values is exceeded, the missing values for that stratum of data are not replaced with imputed values.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT imputed.tsv
# PARAMETER imputation.method [mean, median, knn] DEFAULT knn (Imputation method)
# PARAMETER number.of.neighbors INTEGER FROM 1 TO 100000 DEFAULT 5 (Number of neighbors to use for knn imputation)
# PARAMETER missing.row.max PERCENT DEFAULT 20 (Maximum number of missing values on a row)
# PARAMETER missing.column.max PERCENT DEFAULT 20 (Maximum number of missing values on a column)


# Imputation of missing values by mean or median
# JTT 22.6.2006

# Renaming variables
method<-imputation.method
K<-number.of.neighbors
rmax<-missing.row.max/100
cmax<-missing.column.max/100

# Loads the file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]
annotations<-dat[,grep("annotations", names(dat))]
A<-dat[,grep("average", names(dat))]

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
write.table(dat4, file=("imputed.tsv"), sep="\t", row.names=T, col.names=T, quote=F)