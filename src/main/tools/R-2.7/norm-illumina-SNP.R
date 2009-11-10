# ANALYSIS Normalisation/"Illumina SNP arrays" (Illumina SNP array preprocessing. Input should be a tab-delimited text 
# file with genotype calls. Typically such a file is created using GenCall software from Illumina.) 
# INPUT GENERAL chip.txt OUTPUT normalized.tsv, phenodata.tsv 


# Illumina SNP array normalization
# JTT 22.10.2008

# Reading data
firstfield <- scan(dir(), what = "", sep = ",", flush = TRUE, quiet = TRUE, blank.lines.skip = FALSE, multi.line = FALSE)
skip <- grep("[Data]", firstfield, fixed = TRUE)
samples <- read.table(dir(), skip = skip, header = TRUE, sep = "\t", as.is = TRUE, check.names = FALSE, colClasses = "character")

# Transforming the data to a wide data frame
ids<- unique(samples$"Sample ID")
dat2<-matrix(nrow=nrow(samples)/length(ids), ncol=length(ids), data=NA)
for(i in 1:length(ids)) {
   cursam<-samples[samples$"Sample ID"==ids[i],]
   dat2[,i]<-as.numeric(as.factor(paste(cursam$"Allele1 - Top", cursam$"Allele2 - Top", sep="")))
}
dat2<-data.frame(dat2)
rownames(dat2)<-cursam$"SNP Name"
colnames(dat2)<-paste("chip.", ids, sep="")

# Writes out a phenodata table
chiptype<-"Illumina"
sample<-colnames(dat2)
group<-c(rep("", ncol(dat2)))
training<-c(rep("", ncol(dat2)))
time<-c(rep("", ncol(dat2)))
random<-c(rep("", ncol(dat2)))
write.table(data.frame(sample=sample, chiptype=chiptype, group=group, training=training), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Writing data to disk
write.table(dat2, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
