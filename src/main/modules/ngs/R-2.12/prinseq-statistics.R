# TOOL prinseq-statistics.R: "Reads file statistics" (Calculates general statistics of the reads in the given file. This tool is based on the PRINSEQ program.)
# INPUT fastqfile: "Input reads file" TYPE GENERIC
# OUTPUT reads-stats.tsv 
# OUTPUT OPTIONAL reads-stats.html 
# PARAMETER OPTIONAL input.mode: "Input file format" TYPE [ fq: "FASTQ", fa: "FASTA"] DEFAULT fq (Define the file format of the reads file)

# KM 17.1.2012

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("fastqfile")

#system("
#wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.17.3.tar.gz
#tar zxf prinseq-lite-0.17.3.tar.gz
# ")

# binary
binary.stats <- c(file.path(chipster.tools.path, "prinseq", "prinseq-lite.pl" ))
#binary.stats <- c(file.path("perl prinseq-lite-0.17.3/prinseq-lite.pl"))
# command to generate graph file
system('printf "%s\t%s\t%s\n"  Class Feature Value > reads-stats.tsv')
if (input.mode == "fq") {
	command.stats <- paste("perl", binary.stats, " -fastq fastqfile -out_good null -out_bad null -stats_all >> reads-stats.tsv")
}
if (input.mode == "fa") {
	command.stats <- paste("perl", binary.stats, " -fasta fastqfile -out_good null -out_bad null -stats_all >> reads-stats.tsv")
}

#stop(paste('CHIPSTER-NOTE: ', command.stats))

ret <- system(command.stats)
if (ret > 0) {
	stop('Unsupported input file type, please see tool output for more details.')
}

# commands to generate graph file
if (input.mode == "fq") {
   command.graph <- paste("perl", binary.stats, " -fastq fastqfile -out_good null -out_bad null -graph_data tmp_graph_file ")
}
if (input.mode == "fa") {
   command.graph <- paste("perl", binary.stats, " -fasta fastqfile -out_good null -out_bad null -graph_data tmp_graph_file ")
}
system(command.graph)

## create html file
binary.graph <- c(file.path(chipster.tools.path, "prinseq", "prinseq-graphs.pl"))
command.graph <- paste("perl", binary.graph, " -i tmp_graph_file -html_all -o reads-stats")
system(command.graph)

