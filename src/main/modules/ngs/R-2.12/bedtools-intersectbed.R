# TOOL bedtools-intersectbed.R: "Bedtools intersectBed" (Report overlaps between two feature files.)
# INPUT a TYPE GENERIC
# INPUT b TYPE GENERIC
# OUTPUT result.txt 
# PARAMETER abam: "File A is BAM format" TYPE [yes, no] DEFAULT no (Select yes if file A is BAM format.)
# PARAMETER OPTIONAL ubam: "Write uncompressed BAM output" TYPE [yes, no] DEFAULT no (Write uncompressed BAM output. Default is to write compressed BAM.)
# PARAMETER OPTIONAL bed: "Write output as BED" TYPE [yes, no] DEFAULT no (When using BAM input, write output as BED. The default is to write output in BAM when using BAM input.
# PARAMETER OPTIONAL wa: "Write the original entry in A for each overlap" TYPE [yes, no] DEFAULT no (Write the original entry in A for each overlap.)
# PARAMETER OPTIONAL wb: "Write the original entry in B for each overlap" TYPE [yes, no] DEFAULT no (Write the original entry in B for each overlap.)
# PARAMETER OPTIONAL wo: "Write the original A and B entries plus the number of base pairs of overlap between the two feature" TYPE [yes, no] DEFAULT no
# PARAMETER OPTIONAL wao: "Write the original A and B entries plus the number of base pairs of overlap between the two features" TYPE [yes, no] DEFAULT no
# PARAMETER OPTIONAL u: "Write the original A entry once if any overlaps found in B" TYPE [yes, no] DEFAULT no (Write the original A entry once if any overlaps found in B)
# PARAMETER OPTIONAL c: "For each entry in A, report the number of overlaps with B" TYPE [yes, no] DEFAULT no (For each entry in A, report the number of overlaps with B)
# PARAMETER OPTIONAL v: "Only report those entries in A that have no overlaps with B" TYPE [yes, no] DEFAULT no (Only report those entries in A that have no overlaps with B)
# PARAMETER OPTIONAL f: "Minimum overlap required as a fraction of A" TYPE TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.000000001 (Minimum overlap required as a fraction of A. Default is 1E-9 (effectively 1bp\))
# PARAMETER OPTIONAL r: "Require that the fraction overlap be reciprocal for A and B" TYPE [yes, no] DEFAULT no (Require that the fraction overlap be reciprocal for A and B. In other words, if minumum overlap is 0.90 and this option os selcted, this requires that B overlap 90% of A and A _also_ overlaps 90% of B.)
# PARAMETER OPTIONAL s: "Force strandedness" TYPE [yes, no] DEFAULT no (Only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.)
# PARAMETER OPTIONAL split: "Treat split BAM or BED12 entries as distinct BED intervals" TYPE [yes, no] DEFAULT no (Treat "split" BAM (i.e., having an “N” CIGAR operation\) or BED12 entries as distinct BED intervals.)

# binary
#binary <- "ls > result.txt"

# run
#system(binary)