# ANALYSIS Test/InputOutput (Empty analysis for testing system wide data transfer.)
# INPUT CDNA input.tsv OUTPUT output.tsv

input<-read.table(c("input.tsv"), header=T, sep="\t")
write.table(c(1, 2, 3), file="output.tsv", quote=FALSE, col.names=FALSE, row.names=FALSE)
