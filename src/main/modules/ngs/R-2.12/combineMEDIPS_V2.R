# TOOL combineMEDIPS_V2.R: "Combine several MEDIPS formatted files" (Combines MEDIPS formatted files.)
# INPUT MEDIPS-input{...}.tsv: "MEDIPS files" TYPE GENERIC
# OUTPUT MEDIPS-input.tsv: "A converted BAM file suitable for MEDIPS analysis"

files<-dir(pattern=".tsv")
for(i in 1:length(files)) {
   dat<-read.table(files[i], header=F, sep="\t")
   write.table(dat, "MEDIPS-input.tsv", append=TRUE, quote=F, sep="\t", row.names=FALSE, col.names=FALSE)
}
