# TOOL convert-bam-to-medips.R: "Convert BAM file to MEDIPS input format" (Converts a BAM file to the MEDIPS input format.)
# INPUT bam{...}.bam: "BAM data diles" TYPE GENERIC 
# OUTPUT MEDIPS-input{...}.tsv: "A converted BAM file suitable for MEDIPS analysis" 
# OUTPUT META phenodata.tsv: "Phenodata file needed for MEDIPS analysis" 

# Reads the input files, converts them to the MEDIPS format, and writes the converted files on the disk
library(Rsamtools)
files<-dir(pattern=".bam")
outfiles<-paste("MEDIPS-input", 1:length(files), ".tsv", sep="")
for(i in 1:length(files)) {
   aln<-readBamGappedAlignments(files[i])
   toMEDIPS<-data.frame(rep(rname(aln)@values, rname(aln)@lengths), start(aln), end(aln), rep(strand(aln)@values, strand(aln)@lengths))
   write.table(toMEDIPS, outfiles[i], col.names=F, row.names=F, quote=F, sep="\t")
   rm(aln, toMEDIPS)
   gc()
}

# Writes out a phenodata table
sample<-outfiles
group<-rep("", length(files))
chiptype<-rep("MeDIP", length(files))
write.table(data.frame(sample=sample, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)



