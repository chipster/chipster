# TOOL concatenate.R: "Concatenate two files" (Reads two files and wites the content into one file)
# INPUT input1: "First file" TYPE GENERIC
# INPUT input2: "Second file" TYPE GENERIC
# OUTPUT result.txt


command.full <- paste("cat input1 input2 > result.txt")

system(command.full)




