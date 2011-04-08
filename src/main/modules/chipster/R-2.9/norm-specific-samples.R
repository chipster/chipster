# ANALYSIS Normalisation/"Normalize to specific samples" (Normalizes data to specific samples. 
# The samples to be normalized should be coded with 1 in one column of the phenodata, whereas the samples
# to be used as reference and normalize against should be coded with 0 in the same column.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT normalized2samples.tsv
# PARAMETER column.to.normalize.by METACOLUMN_SEL DEFAULT group (Phenodata column containing the samples to be normalized)
# PARAMETER scale [log, linear] DEFAULT log (Specifies if the data has been log-transformed or is in linear scale.)

# Normalize the data to specific samples
# JTT 5.12.2008
#
# MG, 23.12.2010
# Modified to handle linear scale data as well

# Loads the data file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
dat2<-dat[,grep("chip", names(dat))]

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Extract the data from the phenodata column
extract<-phenodata[,grep(column.to.normalize.by, colnames(phenodata))]

# Sanity checks
if(length(unique(extract))>2) {
   stop("CHIPSTER-NOTE: You have specified more than two groups! You need to define exactly two groups.")
}
if(max(extract>1)) {
   stop("CHIPSTER-NOTE: The groups should be defined with 0s and 1s! You have numbers larger than 1 in the definitions.")
}

# Extracting the samples
dat0<-as.data.frame(dat2[,which(extract==0)])
dat1<-as.data.frame(dat2[,which(extract==1)])

# Normalization values
nv<-(rowSums(dat0)/ncol(dat0))

	
# Normalization
if (scale == "linear") {
	for(i in 1:ncol(dat2)) {
   		dat2[,i]<-dat2[,i]/nv
	}
}
if (scale == "log") {
	for(i in 1:ncol(dat2)) {
		dat2[,i]<-dat2[,i]-nv
	}
}

# Replace the data in the original object with normalized values
datacols<-grep("chip.", colnames(dat))
for(i in 1:length(datacols)) {
   dat[,(datacols[i])]<-dat2[,i]
}

# Write out data
write.table(data.frame(dat), file="normalized2samples.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# EOF
