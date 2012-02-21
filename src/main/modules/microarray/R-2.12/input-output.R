# TOOL input-output.R: InputOutput (Empty analysis for testing system wide data transfer.)
# INPUT input.tsv: input.tsv TYPE GENERIC 
# OUTPUT output.tsv: output.tsv 

#input<-read.table(c("input.tsv"), header=T, sep="\t")
#system("mv input.tsv output.tsv")
write.table(c(1, 2, 3), file="output.tsv", quote=FALSE, col.names=FALSE, row.names=FALSE)
