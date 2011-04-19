# ANALYSIS Test/SimpleInputOutput (Send prefix + input back as output)
# INPUT GENERIC input.dat OUTPUT output.dat

input<-readLines("input.dat")
output<-paste("INPUT WAS:", input)
writeLines(output, "output.dat")

