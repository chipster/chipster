# TOOL prinseq-statistics.R: "FASTQ statistics" (Calculates genral statistics of the reads in the given FASTQ file. This tool is based on the PRINSEQ program.)
# INPUT fastqfile: "Input sequence set" TYPE GENERIC
# OUTPUT reads-stats.tsv 

## OUTPUT reads-stats.html 

# KM 17.1.2012

# check out if the file is compressed and if so unzip it
system("file fastqfile > file_info")
system("grep gzip file_info > is_gzip")
system("[ -s is_gzip ] && mv fastqfile reads.gz ; gzip -d reads.gz ; mv reads fastqfile")

system("
wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.17.3.tar.gz
tar zxf prinseq-lite-0.17.3.tar.gz
 ")

# binary
#binary.stats <- c(file.path(chipster.tools.path, "prinseq", "bin", "prinseq-lite.pl"))
binary.stats <- c(file.path("perl prinseq-lite-0.17.3/prinseq-lite.pl"))
# command to generate graph file
system('printf "%s\t%s\t%s\n"  Class Feature Value > reads-stats.tsv')
command.stats <- paste(binary.stats, " -fastq fastqfile -out_good null -out_bad null -stats_all >> reads-stats.tsv")
# run
#stop(paste('CHIPSTER-NOTE: ', command.stats))
ret <- system(command.stats)
if (ret > 0) {
	stop('Unsupported input file type, please see tool output for more details.')
}

## commands to generate graph file
#command.graph <- paste(binary.stats, " -fastq reads.fastq -out_good null -out_bad null -graph_data tmp_grap_file ")
#system(command.stats)

## create html file
#binary.graph <- c(file.path(chipster.tools.path, "prinseq", "bin", "prinseq-graphs.pl"))
#command.graph <- paste(binary.graph, " -fastq reads.fastq -i tmp_graph_file -html_all -o reads-stats.html")
#system(command.graph)


