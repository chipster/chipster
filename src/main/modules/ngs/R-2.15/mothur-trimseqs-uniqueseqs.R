# TOOL mothur-trimseqs-uniqueseqs.R: "Trim and filter sequences with Mothur" (Removes primers and barcodes, trims and filters reads for several criteria, and removes duplicate reads. This tool is based on the Mothur package.)
# INPUT reads.fasta: "FASTA file" TYPE FASTA
# INPUT reads.oligos: "Oligos" TYPE MOTHUR_OLIGOS
# INPUT OPTIONAL reads.qual: "QUAL file" TYPE GENERIC
# OUTPUT OPTIONAL reads.trim.unique.fasta
# OUTPUT OPTIONAL reads.trim.unique.qual
# OUTPUT OPTIONAL reads.groups
# OUTPUT OPTIONAL reads.trim.names
# OUTPUT OPTIONAL summary.trim.unique.tsv
# PARAMETER OPTIONAL flip: "Use reverse complement" TYPE [yes, no] DEFAULT no (Use reverse complement of the sequences.)
# PARAMETER OPTIONAL qaverage: "Minimum average quality of sequence" TYPE INTEGER FROM 0 TO 40 (Minimum average quality of the sequence. Sequences that have a lower average quality are dropped.)
# PARAMETER OPTIONAL qwindowaverage: "Minimum average quality of window" TYPE INTEGER (Minimum average quality score allowed over a window)
# PARAMETER OPTIONAL qwindowsize: "Window size" TYPE INTEGER (Number of bases in a window)
# PARAMETER OPTIONAL qstepsize: "Window step size" TYPE INTEGER (Number of bases to move the window over.)
# PARAMETER OPTIONAL maxambig: "Maximum ambiguous bases" TYPE INTEGER FROM 0 TO 10 (Maximum number of ambiguous bases allowed in any sequence)
# PARAMETER OPTIONAL maxhomop: "Maximum homopolymer length" TYPE INTEGER FROM 0 TO 50 (Maximum length of a homopolymer allowed in any sequence)
# PARAMETER OPTIONAL minlength: "Minimum sequence length" TYPE INTEGER FROM 0 TO 1000 (Minimum length of an allowed sequence)
# PARAMETER OPTIONAL maxlength: "Maximum sequence length" TYPE INTEGER FROM 0 TO 1000 (Maximum length of an allowed sequence)
# PARAMETER OPTIONAL pdiffs: "Maximum differences to primer sequences" TYPE INTEGER FROM 0 TO 10 (Maximum number of allowed differences to primer sequences)
# PARAMETER OPTIONAL bdiffs: "Maximum differences to barcode sequences" TYPE INTEGER FROM 0 TO 10 (Maximum number of allowed differences to barcode sequences)

# AMS 05.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# Add options
trimseqs.options <- ""
trimseqs.options <- paste(trimseqs.options, "trim.seqs(fasta=reads.fasta, oligos=reads.oligos")
if (file.exists("reads.qual")){
	trimseqs.options <- paste(trimseqs.options, " qfile=reads.qual", sep=",")
}
if (flip == "yes"){
	trimseqs.options <- paste(trimseqs.options, " flip=T", sep=",")
}
if (!is.na(qaverage)){
	trimseqs.options <- paste(trimseqs.options, ", qaverage=", qaverage, sep="")
}
if (!is.na(qwindowaverage)){
	trimseqs.options <- paste(trimseqs.options, ", qwindowaverage=", qwindowaverage, sep="")
}
if (!is.na(qwindowsize)){
	trimseqs.options <- paste(trimseqs.options, ", qwindowsize=", qwindowsize, sep="")
}
if (!is.na(qstepsize)){
	trimseqs.options <- paste(trimseqs.options, ", qstepsize=", qstepsize, sep="")
}
if (!is.na(maxambig)){
	trimseqs.options <- paste(trimseqs.options, ", maxambig=", maxambig, sep="")
}
if (!is.na(maxhomop)){
	trimseqs.options <- paste(trimseqs.options, ", maxhomop=", maxhomop, sep="")
}
if (!is.na(minlength)){
	trimseqs.options <- paste(trimseqs.options, ", minlength=", minlength, sep="")
}
if (!is.na(maxlength)){
	trimseqs.options <- paste(trimseqs.options, ", maxlength=", maxlength, sep="")
}
if (!is.na(pdiffs)){
	trimseqs.options <- paste(trimseqs.options, ", pdiffs=", pdiffs, sep="")
}
if (!is.na(bdiffs)){
	trimseqs.options <- paste(trimseqs.options, ", bdiffs=", bdiffs, sep="")
}
trimseqs.options <- paste(trimseqs.options, ")", sep="")

#stop(paste('CHIPSTER-NOTE: ', trimseqs.options))

# Write batch file
write(trimseqs.options, "trim.mth", append=F)
write("unique.seqs(fasta=reads.trim.fasta)", "trim.mth", append=T)

# command
command <- paste(binary, "trim.mth")

# run
system(command)


# batch file
write("summary.seqs(fasta=reads.trim.unique.fasta)", "summary.mth", append=F)

# command
command <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command)

# Make reads.trim.unique.qual
if (file.exists("reads.trim.qual")){
	system("grep '>' reads.trim.unique.fasta | cut -c 2- > reads.trim.unique.list")
	system("perl -ne 'if(/^>(\\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' reads.trim.unique.list reads.trim.qual > reads.trim.unique.qual")
}

# Postprocess output files
system("grep -A 9 Start log_raw.txt > summary.trim.unique.tsv")
#system("mv reads.trim.names reads.trim.names.txt")
#system("mv reads.groups reads.groups.txt")








