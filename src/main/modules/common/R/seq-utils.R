# Adds "chr" to the beginning of each line of a fasta file that starts with a number or with X, Y, Z, W or M
#
addChrToFasta <- function(input, output){
	system(paste("sed /^.[0-9,X,Y,Z,W,M]/s/\">\"/\">chr\"/  <", input, ">", output))
}