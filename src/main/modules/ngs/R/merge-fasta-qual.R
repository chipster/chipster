# TOOL merge-fasta-qual.R: "Merge FASTA or QUAL files" (Merge two or more FASTA or QUAL files together. You should take care to only merge files of the same type.)
# INPUT file{...}.tmp: "Files to be merged" TYPE GENERIC
# OUTPUT OPTIONAL merged.fasta
# OUTPUT OPTIONAL merged.qual
# PARAMETER type: "File type" TYPE [FASTA, QUAL] DEFAULT FASTA (Type of files to be merged.)

# AMS 19.06.2013

if (type == "FASTA"){
	system("cat *.tmp > merged.fasta")
}
if (type == "QUAL"){
	system("cat *.tmp > merged.qual")
}
