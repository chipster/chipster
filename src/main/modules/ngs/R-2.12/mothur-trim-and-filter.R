# TOOL mothur-trim-and-filter.R: "Trim and filter reads with Mothur" (Trim and filter Reads. This tool is based on the Mothur package.)
# INPUT reads.fasta: "FASTA file" TYPE GENERIC
# INPUT reads.oligos: "Oligos" TYPE GENERIC
# INPUT OPTIONAL reads.qual: "QUAL file" TYPE GENERIC
# OUTPUT OPTIONAL reads.trim.unique.fasta
# OUTPUT OPTIONAL reads.trim.unique.qual
# OUTPUT OPTIONAL reads.trim.qual
# OUTPUT OPTIONAL reads.groups
# OUTPUT OPTIONAL reads.trim.names
# OUTPUT OPTIONAL summary.trim.unique.tsv
# PARAMETER OPTIONAL flip: "Use reverse complement" TYPE [yes, no] DEFAULT no (Use reverse complement of the sequences.)
# PARAMETER OPTIONAL qaverage: "Minimum average quality of sequence" TYPE INTEGER FROM 0 TO 40 DEFAULT 20 (Minimum average quality of the sequence. Sequences that have a lower average quality are dropped.)
# PARAMETER OPTIONAL qwindowaverage: "Minimum average quality of Window" TYPE INTEGER DEFAULT 20 (Minimum average quality score allowed over a window)
# PARAMETER OPTIONAL qwindowsize: "Window size" TYPE INTEGER DEFAULT 50 (Number of bases in a window)
# PARAMETER OPTIONAL qstepsize: "Window step size" TYPE INTEGER DEFAULT 1 (Number of bases to move the window over.)
# PARAMETER OPTIONAL maxambig: "Maximum ambiguos bases" TYPE INTEGER FROM 0 TO 10 DEFAULT 0 (Maximum number of ambiguous bases allowed in any sequence)
# PARAMETER OPTIONAL maxhomop: "Maximum homopolymer length" TYPE INTEGER FROM 0 TO 50 DEFAULT 8 (Maximum length of a homopolymere allowed in any sequence)
# PARAMETER OPTIONAL minlength: "Minimum sequence length" TYPE INTEGER FROM 0 TO 1000 DEFAULT 250 (Minimum length of an allowed sequence)
# PARAMETER OPTIONAL maxlength: "Maximum sequence length" TYPE INTEGER FROM 0 TO 1000 DEFAULT 350 (Maximum length of an allowed sequence)
# PARAMETER OPTIONAL pdiffs: "Maximum differences to primer sequences" TYPE INTEGER FROM 0 TO 10 DEFAULT 0 (Maximum number of allowed differences to primer sequences)
# PARAMETER OPTIONAL bdiffs: "Maximum differences to barcode sequences" TYPE INTEGER FROM 0 TO 10 DEFAULT 0 (Maximum number of allowed differences to barcode sequences)

# AMS 05.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "1.28.0", "mothur"))

# Add options
trimseqs.options <- ""
trimseqs.options <- paste(trimseqs.options, "trim.seqs(fasta=reads.fasta, oligos=reads.oligos")
if (file.exists("reads.qual")){
	trimseqs.options <- paste(trimseqs.options, " qfile=reads.qual", sep=",")
}
if (flip == "yes"){
	trimseqs.options <- paste(trimseqs.options, " flip=T", sep=",")
}
trimseqs.options <- paste(trimseqs.options, ", qaverage=", qaverage, sep="")
trimseqs.options <- paste(trimseqs.options, ", qwindowaverage=", qwindowaverage, sep="")
trimseqs.options <- paste(trimseqs.options, ", qwindowsize=", qwindowsize, sep="")
trimseqs.options <- paste(trimseqs.options, ", qstepsize=", qstepsize, sep="")
trimseqs.options <- paste(trimseqs.options, ", maxambig=", maxambig, sep="")
trimseqs.options <- paste(trimseqs.options, ", maxhomop=", maxhomop, sep="")
trimseqs.options <- paste(trimseqs.options, ", minlength=", minlength, sep="")
trimseqs.options <- paste(trimseqs.options, ", maxlength=", maxlength, sep="")
trimseqs.options <- paste(trimseqs.options, ", pdiffs=", pdiffs, sep="")
trimseqs.options <- paste(trimseqs.options, ", bdiffs=", bdiffs, sep="")
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

# Post process output
#system("mv reads.summary summary.tsv")
system("grep -A 9 Start log_raw.txt > summary.trim.unique.tsv")

# Post process output

# Make reads.trim.unique.qual
if (file.exists("reads.trim.qual")){
	system("grep '>' reads.trim.unique.fasta | cut -c 2- > reads.trim.unique.list")
	system("perl -ne 'if(/^>(\\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' reads.trim.unique.list reads.trim.qual > reads.trim.unique.qual")
}
#system("mv reads.trim.unique.fasta trimmed.unique.fasta")
#system("mv reads.trim.unique.fasta trimmed.unique.fasta")
#system("grep -A 9 Start log_raw.txt > log.txt")




