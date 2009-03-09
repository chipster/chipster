# Download promoter sequences from UCSC :
# http://hgdownload.cse.ucsc.edu/downloads.html
# and from:
# http://www.ncbi.nlm.nih.gov/CBBresearch/Landsman/Cell_cycle_data/upstream_seq.html

# Loads the libraries
library(Biostrings)

# Human
fas<-readFASTA("upstream1000.fa")
v<-rep(NA, length(fas))
for(i in 1:length(fas)) {
   v[i]<-fas[[i]]$desc
}
v2<-substr(v, 2, 12)
v3<-gsub(x=v2, pattern="_[a-z]", "")
s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
for(i in 1:length(fas)) {
   s[i,]<-fas[[i]]$seq
}
names(s)<-"Sequence"
write.table(data.frame(RefSeq=v3, Sequence=s), col.names=T, row.names=F, sep="\t", quote=F, file="UCSC_Build_36.1_hg18_upstream1000.tsv")
rm(list=objects())

fas<-readFASTA("upstream2000.fa")
v<-rep(NA, length(fas))
for(i in 1:length(fas)) {
   v[i]<-fas[[i]]$desc
}
v2<-substr(v, 2, 12)
v3<-gsub(x=v2, pattern="_[a-z]", "")
s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
for(i in 1:length(fas)) {
   s[i,]<-fas[[i]]$seq
}
names(s)<-"Sequence"
write.table(data.frame(RefSeq=v3, Sequence=s), col.names=T, row.names=F, sep="\t", quote=F, file="UCSC_Build_36.1_hg18_upstream2000.tsv")
rm(list=objects())

fas<-readFASTA("upstream5000.fa")
v<-rep(NA, length(fas))
for(i in 1:length(fas)) {
   v[i]<-fas[[i]]$desc
}
v2<-substr(v, 2, 12)
v3<-gsub(x=v2, pattern="_[a-z]", "")
s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
for(i in 1:length(fas)) {
   s[i,]<-fas[[i]]$seq
}
names(s)<-"Sequence"
write.table(data.frame(RefSeq=v3, Sequence=s), col.names=T, row.names=F, sep="\t", quote=F, file="UCSC_Build_36.1_hg18_upstream5000.tsv")
rm(list=objects())


# Mouse
fas<-readFASTA("upstream1000.fa")
v<-rep(NA, length(fas))
for(i in 1:length(fas)) {
   v[i]<-fas[[i]]$desc
}
v2<-substr(v, 2, 12)
v3<-gsub(x=v2, pattern="_[a-z]", "")
s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
for(i in 1:length(fas)) {
   s[i,]<-fas[[i]]$seq
}
names(s)<-"Sequence"
write.table(data.frame(RefSeq=v3, Sequence=s), col.names=T, row.names=F, sep="\t", quote=F, file="UCSC_Build_36_mm8_upstream1000.tsv")
rm(list=objects())

fas<-readFASTA("upstream2000.fa")
v<-rep(NA, length(fas))
for(i in 1:length(fas)) {
   v[i]<-fas[[i]]$desc
}
v2<-substr(v, 2, 12)
v3<-gsub(x=v2, pattern="_[a-z]", "")
s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
for(i in 1:length(fas)) {
   s[i,]<-fas[[i]]$seq
}
names(s)<-"Sequence"
write.table(data.frame(RefSeq=v3, Sequence=s), col.names=T, row.names=F, sep="\t", quote=F, file="UCSC_Build_36_mm8_upstream2000.tsv")
rm(list=objects())

fas<-readFASTA("upstream5000.fa")
v<-rep(NA, length(fas))
for(i in 1:length(fas)) {
   v[i]<-fas[[i]]$desc
}
v2<-substr(v, 2, 12)
v3<-gsub(x=v2, pattern="_[a-z]", "")
s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
for(i in 1:length(fas)) {
   s[i,]<-fas[[i]]$seq
}
names(s)<-"Sequence"
write.table(data.frame(RefSeq=v3, Sequence=s), col.names=T, row.names=F, sep="\t", quote=F, file="UCSC_Build_36_mm8_upstream5000.tsv")
rm(list=objects())


# Rat
fas<-readFASTA("upstream1000.fa")
v<-rep(NA, length(fas))
for(i in 1:length(fas)) {
   v[i]<-fas[[i]]$desc
}
v2<-substr(v, 2, 12)
v3<-gsub(x=v2, pattern="_[a-z]", "")
s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
for(i in 1:length(fas)) {
   s[i,]<-fas[[i]]$seq
}
names(s)<-"Sequence"
write.table(data.frame(RefSeq=v3, Sequence=s), col.names=T, row.names=F, sep="\t", quote=F, file="UCSC_rn4_upstream1000.tsv")
rm(list=objects())

fas<-readFASTA("upstream2000.fa")
v<-rep(NA, length(fas))
for(i in 1:length(fas)) {
   v[i]<-fas[[i]]$desc
}
v2<-substr(v, 2, 12)
v3<-gsub(x=v2, pattern="_[a-z]", "")
s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
for(i in 1:length(fas)) {
   s[i,]<-fas[[i]]$seq
}
names(s)<-"Sequence"
write.table(data.frame(RefSeq=v3, Sequence=s), col.names=T, row.names=F, sep="\t", quote=F, file="UCSC_rn4_upstream2000.tsv")
rm(list=objects())

fas<-readFASTA("upstream5000.fa")
v<-rep(NA, length(fas))
for(i in 1:length(fas)) {
   v[i]<-fas[[i]]$desc
}
v2<-substr(v, 2, 12)
v3<-gsub(x=v2, pattern="_[a-z]", "")
s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
for(i in 1:length(fas)) {
   s[i,]<-fas[[i]]$seq
}
names(s)<-"Sequence"
write.table(data.frame(RefSeq=v3, Sequence=s), col.names=T, row.names=F, sep="\t", quote=F, file="UCSC_rn4_upstream5000.tsv")
rm(list=objects())


# Yeast
fas<-readFASTA("yeast500.fas")
v<-rep(NA, length(fas))
for(i in 1:length(fas)) {
   v[i]<-fas[[i]]$desc
}
v2<-substr(v, 2, 12)
v3<-gsub(x=v2, pattern="_[a-z]", "")
s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
for(i in 1:length(fas)) {
   s[i,]<-fas[[i]]$seq
}
names(s)<-"Sequence"
write.table(data.frame(RefSeq=v3, Sequence=s), col.names=T, row.names=F, sep="\t", quote=F, file="NCBI_sc_upstream500.tsv")
rm(list=objects())

fas<-readFASTA("yeast1000.fas")
v<-rep(NA, length(fas))
for(i in 1:length(fas)) {
   v[i]<-fas[[i]]$desc
}
v2<-substr(v, 2, 12)
v3<-gsub(x=v2, pattern="_[a-z]", "")
s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
for(i in 1:length(fas)) {
   s[i,]<-fas[[i]]$seq
}
names(s)<-"Sequence"
write.table(data.frame(RefSeq=v3, Sequence=s), col.names=T, row.names=F, sep="\t", quote=F, file="NCBI_sc_upstream1000.tsv")
rm(list=objects())

fas<-readFASTA("yeast2500.fas")
v<-rep(NA, length(fas))
for(i in 1:length(fas)) {
   v[i]<-fas[[i]]$desc
}
v2<-substr(v, 2, 12)
v3<-gsub(x=v2, pattern="_[a-z]", "")
s<-data.frame(matrix(nrow=length(fas), ncol=1, data=NA))
for(i in 1:length(fas)) {
   s[i,]<-fas[[i]]$seq
}
names(s)<-"Sequence"
write.table(data.frame(RefSeq=v3, Sequence=s), col.names=T, row.names=F, sep="\t", quote=F, file="NCBI_sc_upstream2500.tsv")
rm(list=objects())
