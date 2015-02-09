# Adds "chr" to the each chromosome name starting with number or X or Y, Changes MT to chrM
# It uses samtools reheader command
#
addChrToBAM <- function(input, output){
	samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
		system(paste(samtools.binary, "view -H", input, "| sed -e 's/SN:\\([0-9XY]\\)/SN:chr\\1/' -e 's/SN:MT/SN:chrM/' |", samtools.binary, "reheader -", input, ">", output))
}							
