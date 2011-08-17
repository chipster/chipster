# TOOL bedtools-mergebed.R: "BEDTools mergeBed" (mergeBed combines overlapping or book-ended (that is, one base pair away\) features in a feature file into a single feature which spans all of the combined features.)
# INPUT file.a: "Input file" TYPE GENERIC
# OUTPUT mergebed.bed
# PARAMETER OPTIONAL s: "Force strandedness" TYPE [yes, no] DEFAULT no (Force strandedness.  That is, only merge features that are the same strand. By default, merging is done without respect to strand.)
# PARAMETER OPTIONAL n: "Report the number of BED entries that were merged" TYPE [yes, no] DEFAULT no (Report the number of BED entries that were merged. Note: 1 is reported if no merging occurred.)
# PARAMETER OPTIONAL d: "Maximum distance between features allowed" TYPE INTEGER DEFAULT 0 (Maximum distance between features allowed for features to be merged. Default is 0. That is, overlapping & book-ended features are merged.)
# PARAMETER OPTIONAL nms: "Report the names of the merged features separated by semicolons" TYPE [yes, no] DEFAULT no (Report the names of the merged features separated by semicolons.)
# PARAMETER OPTIONAL score: "Report the scores of the merged features" TYPE [yes, no] DEFAULT no (Report the scores of the merged features. Specify one of the following options for reporting scores: sum, min, max, mean, median, mode, antimode, collapse (i.e., print a semicolon-separated list\).)
# PARAMETER OPTIONAL score.type: "Score report type" TYPE [sum, min, max, mean, median, mode, antimode, collapse] DEFAULT sum ()

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "mergeBed"))

# options
options <- paste("")
if (s == "yes") {options <- paste(options, "-s")}
if (n == "yes") {options <- paste(options, "-n")}
options <- paste(options, "-d", d)
if (nms == "yes") {options <- paste(options, "-nms")}
if (score == "yes") {options <- paste(options, "-score", score.type)}

# input files
options <- paste(options, "-i file.a")

# command
command <- paste(binary, options, "> mergebed.bed")

# run
system(command)