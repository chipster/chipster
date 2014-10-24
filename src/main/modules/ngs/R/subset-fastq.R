# TOOL subset-fastq.R: "Make a subset of FASTQ" (Returns a subsample of N reads from a FASTQ file. Tool is based on the seqtk package.)
# INPUT reads.fastq: "FASTQ file" TYPE GENERIC
# OUTPUT subset.fastq.gz 
# PARAMETER n.seq: "Size of subset" TYPE INTEGER DEFAULT 100000 (Number of reads to return from the FASTQ file.)
# PARAMETER seed: "Random seed" TYPE INTEGER DEFAULT 11 (Random seed for the sampling. When using paired-end data, use same random seed to keep pairing.)

# AMS 14.5.2012
# EK 15.5.2012 added unzipping
# AMS 24.9.2014 added zipping the result file
# AMS 13.10.2014 changed to use seqtk

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")

# binary
seqtk.binary <- file.path(chipster.tools.path, "seqtk", "seqtk")

seed.option <- paste ("-s",as.character(seed), sep="")

# command
command <- paste(seqtk.binary, "sample", seed.option, "reads.fastq", n.seq, "> subset.fastq")

# run
#stop(paste('CHIPSTER-NOTE: ', command))
system(command)
system("gzip subset.fastq")