# TOOL mothur-sffinfo.R: "Extract sequences from a SFF file" (Extract a FASTA sequence file and a QUAL quality file from a SFF file. Sequences can optionally be trimmed for quality. This tool is based on the Mothur package.)
# INPUT reads.sff: "SFF file" TYPE GENERIC
# OUTPUT OPTIONAL reads.fasta
# OUTPUT OPTIONAL reads.qual
# OUTPUT OPTIONAL reads.raw.fasta
# OUTPUT OPTIONAL reads.raw.qual
# PARAMETER OPTIONAL trim: "Trim reads for quality" TYPE [yes, no] DEFAULT no (Trim sequences and quality scores to the clipQualLeft and clipQualRight values.)

# AMS 19.06.2013

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# Options
sffinfo.options <- ""
sffinfo.options <- paste(sffinfo.options, "sffinfo(sff=reads.sff")
if (trim == "no"){
	sffinfo.options <- paste(sffinfo.options, " trim=F", sep=",")
}
sffinfo.options <- paste(sffinfo.options, ")", sep="")

#stop(paste('CHIPSTER-NOTE: ', ssfinfo.options))

# Write batch file
write(sffinfo.options, "sffinfo.mth", append=F)

# command
command <- paste(binary, "sffinfo.mth")

# run
system(command)

