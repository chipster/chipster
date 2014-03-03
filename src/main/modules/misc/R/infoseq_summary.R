# TOOL infoseq_summary.R: "Sequence file summary" (Tool to calculate basic properties of a sequence file.)
# INPUT input.txt: "Query sequences" TYPE GENERIC
# OUTPUT summary_file.txt

# K.M 28.10.2013

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("input.txt")

# pb settings
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
#infoseq.binary <- file.path(chipster.tools.path, "blast", "/ncbi-blast-2.2.28+", "bin", "infoseq_summary.bash")
infoseq.binary <- file.path(chipster.module.path ,"/shell/infoseq_summary.sh")
command.full <- paste(infoseq.binary, emboss.path, "input.txt > summary_file.txt 2>&1" )
system(command.full)
