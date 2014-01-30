# TOOL simple-input-output.R: SimpleInputOutput (Send prefix + input back as output)
# INPUT input.dat: input.dat TYPE GENERIC 
# OUTPUT output.dat: output.dat 

input<-readLines("input.dat")
output<-paste("INPUT WAS:", input)
writeLines(output, "output.dat")

