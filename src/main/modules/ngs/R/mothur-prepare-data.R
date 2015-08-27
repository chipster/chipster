# TOOL mothur-prepare-data.R: "Mothur - prepare data" (Preprocesses the sequence and the quality data for further analysis in mothur. Select all fasta and quality files you want to process at the same time as input files. This tool appends all sequence files together and does the same for all quality files, also.)
# INPUT file{...}.txt: "Sequence or data quality files" TYPE GENERIC 
# OUTPUT all.fasta
# OUTPUT all.qual
# OUTPUT META phenodata.tsv
# OUTPUT OPTIONAL prepare_data_log.txt
# PARAMETER number.to.sample: "Number to sample" TYPE INTEGER FROM 0 TO 1000000 DEFAULT 0 (Number of sequences to sample from each each sequence file. If the value is left at default, all sequences are selected.)

# JTT 2012-11-01


# Assessory function
f<-function(x) {
	seq(x[2], x[3])
}

# Checking the file types
files<-dir()
type<-rep(NA, length(files))
for(i in 1:length(files)) {
	d<-readLines(files[i], n=2)
	if(!any(names(table(strsplit(d[2], ""))) %in% c("A", "C", "G", "T"))) {
		type[i]<-"qual"
	}
	if(d[1]==".sff") {
		type[i]<-"sff"
	}
	if(any(names(table(strsplit(d[2], ""))) %in% c("A", "C", "G", "T"))) {
		type[i]<-"seq"
	}
}
qual<-files[type=="qual"]
fasta<-files[type=="seq"]

# Process the files, if no sampling is used
if(number.to.sample==0) {
	for(i in 1:length(qual)) {
		d<-readLines(qual[i])
		write(d, "all.qual", append=T)
	}
	
	for(i in 1:length(fasta)) {
		d<-readLines(fasta[i])
		write(d, "all.fasta", append=T)
	}
}

# Process the files, if sampling is used
if(number.to.sample!=0) {
	for(i in 1:length(fasta)) {
		# Sequence files
		d<-readLines(fasta[i])
		s=grep(">", d)
		h=d[s]
		e=c((s-1)[-1], length(d))
		h<-data.frame(h, s, e)
		if(nrow(h)<number.to.sample) {
			stop(paste("There are less sequences than you specified for sampling! Sample ", fasta[i], " contains only ", nrow(h), " sequences. Please use rerun the tool using a smaller number of sequences.",sep=""))
		}
		s<-sample(1:nrow(h), number.to.sample, replace=FALSE)
		hs<-h[s,]
		l<-sort(as.vector(unlist(apply(hs, 1, f))))
		write(d[l], "all.fasta", append=T)
		# Quality files
		qu<-readLines(qual[i])
		write(qu[l], "all.qual", append=T)
	}
}

# Write out phenodata table
sample<-c(" "," ", fasta)
group<-c(rep("", length(sample))) 
oligo<-c("forward", "reverse", rep("barcode", length(fasta))) 
sequence<-rep("", length(sample))
write.table(data.frame(sample=sample, oligo=oligo, sequence=sequence, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)
