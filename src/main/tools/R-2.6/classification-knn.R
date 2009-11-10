# ANALYSIS Clustering/"KNN classification" (K-nearest neighbor classification. If you have a separate test data set,
# you can validate the prediction with it by setting the validation type to predict. This function does not
# perform any gene selection, and the analysis is run for the selected data set.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT knn-cross-validation.tsv
# PARAMETER number.of.nearest.neighbors INTEGER FROM 1 TO 1000 DEFAULT 2 (Number of nearest neighbors)
# PARAMETER number.of.votes INTEGER FROM 0 TO 1000 DEFAULT 2 (Number of votes needed to get a definite answer)
# PARAMETER validation.type [crossvalidate, predict] DEFAULT crossvalidate (Type of analysis)


# KNN classification
# JTT 26.6.2006

# Renaming the variables
knn.type<-validation.type
k.no<-number.of.nearest.neighbors
k.vote<-number.of.votes

# Load the libraries
library(class)

# Loads the data
file<-c("normalized.tsv")
dat<-read.table(file, sep="\t", header=T, row.names=1)

# Reads the phenodata table
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Checks whether the training variable is present if phenodata
# If training part of phenodata has not been filled, but the column (header) is present, 
# all the chips belong to the training set
if(grep("training", names(phenodata))>0) {
   tr<-phenodata$training
} 
if(tr[1]==" " | tr[1]=="" | any(is.na(tr))==T) {
   tr<-rep(1, nrow(phenodata))
}

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Are the parameter values sensical?
if(k.no>length(dat2)) {
   stop("The number of neighbors is larger than the number of chips!")
}
if(k.vote>k.no) {
   stop("The number of votes needed to give a definitive answer is larger than the number of neighbors!")
}

# Which parts of the data are training and test sets?
dat3<-split(as.data.frame(t(dat2)), tr) 
train<-dat3$'1'
test<-dat3$'2'

# Defines the true classification of the training set
cl<-phenodata$group

# Runs the KNN analysis and reports the results
if(knn.type=="crossvalidate") {
   knn.cross<-knn.cv(train=train, cl=cl, k=k.no, l=k.vote)
   # Writes a table of known versus predicted classes
   write.table(data.frame(sample=names(dat2), known.classes=cl, prediction=knn.cross), file="knn-cross-validation.tsv", sep="\t", row.names=F, col.names=T, quote=F)
}

if(knn.type=="predict") {
   cl<-split(cl, tr)$'1'
   knn.predict<-knn(train=train, test=test, cl=cl, k=k.no, l=k.vote)
   # Writes a table of known versus predicted classes 
   write.table(data.frame(sample=names(dat2), known.classes=cl, prediction=knn.cross), file="knn-cross-validation.tsv", sep="\t", row.names=F, col.names=T, quote=F)
}