# TOOL bedtools-bamtobed.R: "BEDTools bamToBed" (Converts BAM alignments to BED6, BED12 or BEDPE format.)
# INPUT file.a: "BAM file" TYPE GENERIC
# OUTPUT bamtobed.bam
# PARAMETER output.type: "Output type" TYPE [BED6, BED12, BEDPE] DEFAULT "BED6" (Select the output type (BED6, BED12 or BEDPE\).)
# PARAMETER OPTIONAL split: "Report split BAM alignments as separate BED entries" TYPE [yes,no] DEFAULT no (Report each portion of a split BAM (i.e., having an N CIGAR operation\) alignment as a distinct BED intervals.)
# PARAMETER OPTIONAL ed: "Use BAM edit distance for BED score" TYPE [yes,no] DEFAULT no (Use BAM edit distance (NM tag\) for BED score. Default for BED is to use mapping quality. Default for BEDPE is to use the minimum of the two mapping qualities for the pair.)
# PARAMETER OPTIONAL color: "Color string" TYPE STRING DEFAULT "255,0,0" (An R,G,B string for the color used with BED12 format. Default is (255,0,0\).)
# PARAMETER OPTIONAL cigar: "Add the CIGAR string" TYPE [yes,no] DEFAULT no (Add the CIGAR string to the BED entry as a 7th column.)

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "bamToBed"))



# run
system(binary)