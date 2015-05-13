# TOOL mothur-combine-results.R: mothur-combine-results (Combine results from two or more sequence processing runs.)
# INPUT file.txt: "Groups and taxonomy files from mothur runs" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata files for the mothur runs" TYPE GENERIC  
# OUTPUT all.tax: all.tax
# OUTPUT all.grp: all.grp
# OUTPUT META phenodata-merged.tsv: phenodata.tsv


# JTT 2012-11-05


# Check the file types
files<-dir()
files<-files[-c(grep("phenodata", files))]
types<-rep(NA, length(files))
for(i in 1:length(files)) {
	d<-readLines(files[i], n=1)
	if(length(grep("Bacteria", d))==0) {
		types[i]<-"tax"
	} else {
		types[i]<-"grp"
	}
}

# Append data files
for(i in 1:length(files)) {
	d<-readLines(files[i])   
	if(types[i]=="tax") {
		write(d, "all.tax", append=TRUE)
	}
	if(types[i]=="grp") {
		write(d, "all.grp", append=TRUE)
	}
}

# Combine phenodata tables
files<-dir(pattern="phenodata")
p<-list()
for(i in 1:length(files)) {
	p[[i]]<-read.table(files[i], header=T, sep="\t")[-c(1,2),]
}
pdata<-NULL 
for(i in 1:length(p)) {
	pdata<-rbind(pdata, p[[i]])
}
write.table(pdata, file="phenodata-merged.tsv", quote=FALSE, sep="\t", na="", row.names=FALSE, col.names=TRUE)