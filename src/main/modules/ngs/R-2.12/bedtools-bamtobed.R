# TOOL bedtools-bamtobed.R: "BEDTools bamToBed" (Converts BAM alignments to BED6, BED12 or BEDPE format.)
# INPUT file.a: "BAM file" TYPE GENERIC
# OUTPUT bamtobed.bam
# PARAMETER OPTIONAL bedpe: "Write BEDPE format" TYPE [yes,no] DEAFAULT no (Write BEDPE format. Requires BAM to be grouped or sorted by query.)
# PARAMETER OPTIONAL bed12: "Write BED12" TYPE [yes,no] DEAFAULT no (Write blocked BED format (aka BED12\).)
# PARAMETER OPTIONAL split: "Report split BAM alignments as separate BED entries" TYPE [yes,no] DEAFAULT no (Report each portion of a split BAM (i.e., having an N CIGAR operation\) alignment as a distinct BED intervals.)
# PARAMETER OPTIONAL ed: "Use BAM edit distance (NM tag) for BED score" TYPE [yes,no] DEAFAULT no (Use BAM edit distance (NM tag) for BED score. Default for BED is to use mapping quality. Default for BEDPE is to use the minimum of the two mapping qualities for the pair.)
# PARAMETER OPTIONAL tag: "Use other NUMERIC BAM alignment tag for BED score" TYPE [yes,no] DEAFAULT no (Use other NUMERIC BAM alignment tag for BED score. Default for BED is to use mapping quality. Disallowed with BEDPE output.)
# PARAMETER OPTIONAL cigar: "Add the CIGAR string" TYPE [yes,no] DEAFAULT no (Add the CIGAR string to the BED entry as a 7th column.)

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "bamToBed"))



# run
system(binary)