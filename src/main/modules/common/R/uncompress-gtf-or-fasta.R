# TOOL uncompress-gtf-or-fasta.R: "Uncompress a file" (Uncompress a GTF or a FASTA file. Most Chipster tools accept compressed files directly, but e.g. Genome Browser requires uncompressed files.)
# INPUT input: "Compressed file" TYPE GENERIC
# OUTPUT OPTIONAL uncompressed.gtf
# OUTPUT OPTIONAL uncompressed.fa
# PARAMETER type: "File type" TYPE [GTF, FASTA] DEFAULT GTF (File type.)

source(file.path(chipster.common.path, "zip-utils.R"))
	unzipIfGZipFile("input")

if (type == "GTF"){
	system("mv input uncompressed.gtf")
	
}else{
	system("mv input uncompressed.fa")
}






