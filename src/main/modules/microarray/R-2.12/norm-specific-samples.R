# TOOL norm-specific-samples.R: "Normalize to specific samples" (Normalizes data to specific samples. The samples to be normalized are coded with 1 in one column of the phenodata. The samples to be normalized to are coded with 0 in the same column.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT normalized2samples.tsv: normalized2samples.tsv 
# PARAMETER column.to.normalize.by: column.to.normalize.by TYPE METACOLUMN_SEL DEFAULT group (Phenodata column containing the samples to be normalized)


# Normalize the data to specific samples
# JTT 5.12.2008

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
   stop("You have specified more than two groups! You need to define exactly two groups.")
}
if(max(extract>1)) {
   stop("The groups should be defined with 0s and 1s! You have numbers larger than 1 in the definitions.")
}

# Extracting the samples
dat0<-as.data.frame(dat2[,which(extract==0)])
dat1<-as.data.frame(dat2[,which(extract==1)])

# Normalization values
nv<-(rowSums(dat0)/ncol(dat0))

# Normalization
for(i in 1:ncol(dat2)) {
   dat2[,i]<-dat2[,i]-nv
}

# Replace the data in the original object with normalized values
datacols<-grep("chip.", colnames(dat))
for(i in 1:length(datacols)) {
   dat[,(datacols[i])]<-dat2[,i]
}

# Write out data
write.table(data.frame(dat), file="normalized2samples.tsv", sep="\t", row.names=T, col.names=T, quote=F)
