# TOOL picard-collectmultiplemetrics.R: "Collect multiple metrics from BAM" (Takes an input BAM and runs some Picard metrics modules. This tool is based on the Picard Tools package.)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# OUTPUT OPTIONAL alignment_summary_metrics.tsv           
# OUTPUT OPTIONAL insert_size_metrics.tsv
# OUTPUT OPTIONAL reports.pdf 

# 2015.09.09 AMS 

picard.binary <- file.path(chipster.tools.path, "picard-tools", "picard.jar")
rscript.path <- file.path(chipster.tools.path, "R-3.0.2", "bin")


# run
set.path <-paste(sep="", "PATH=", rscript.path, ":$PATH")
picard.command <- paste("java -Xmx2g -jar", picard.binary, "CollectMultipleMetrics INPUT=alignment.bam OUTPUT=cmm")
command <- paste("bash -c '", set.path, picard.command, "'")
system(command)

# alignment_summary_metrics
system("grep -A2 CATEGORY cmm.alignment_summary_metrics > alignment_summary_metrics.tsv")

# insert_size_metrics
system("grep -A2 CATEGORY cmm.insert_size_metric > insert_size_metric.tsv")

# Join the PDFs 
system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=reports.pdf *.pdf")

# Handle output names
#source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
#inputnames <- read_input_definitions()

# Make a matrix of output names
#outputnames <- matrix(NA, nrow=1, ncol=2)
#outputnames[1,] <- c("unique_alignments.bam", paste(inputnames$alignment.bam))

# Write output definitions file
#write_output_definitions(outputnames)
