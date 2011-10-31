# TOOL bedtools-genomecoveragebed.R: "Genome coverage BED" (Compute the coverage of a feature file among a genome. This tool is based on the BEDTools package.)
# INPUT file.a: "Input file" TYPE GENERIC
# INPUT file.b: "Genome file" TYPE GENERIC
# OUTPUT genomecoveragebed.txt 
# PARAMETER ibam: "Input file is BAM format" TYPE [yes, no] DEFAULT no (The input file is in BAM format. Note: BAM must be sorted by position.)
# PARAMETER output: "Output type" TYPE [histogram, depth, BedGraph, BedGraph-all] DEFAULT histogram (Set the output type. Depth = the depth at each genome position. BedGraph = depth in BedGraph format. BedGraph-all = BedGraph format with regions with zero coverage also reported.)
# PARAMETER OPTIONAL split: "Treat split BAM or BED12 entries as distinct BED intervals" TYPE [yes, no] DEFAULT no (Treat split BAM or BED12 entries as distinct BED intervals when computing coverage. For BAM files, this uses the CIGAR N and D operations to infer the blocks for computing coverage. For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds fields (i.e., columns 10,11,12\).)
# PARAMETER OPTIONAL strand: "Calculate coverage from specific strand" TYPE [both, +, -] DEFAULT both (Calculate coverage of intervals from a specific strand. With BED files, requires at least 6 columns (strand is column 6\).)
# PARAMETER OPTIONAL max: "Combine positions by depth" TYPE [yes, no] DEFAULT no (Combine all positions with a depth >= max depth into a single bin in the histogram. Option is relevant only for histogram output.)
# PARAMETER OPTIONAL maxdepth: "Max depth" TYPE INTEGER DEFAULT 1 ()

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "genomeCoverageBed"))

# optional options
options <- paste("")
if (output == "depth") {options <- paste(options, "-d")}
if (output == "BedGraph") {options <- paste(options, "-bg")}
if (output == "BedGraph-all") {options <- paste(options, "-bga")}
if (split == "yes") {options <- paste(options, "-split")}
if (strand == "+") {options <- paste(options, "-strand +")}
if (strand == "-") {options <- paste(options, "-strand -")}
if (max == "yes" && output == "histogram"){options <- paste(options, "-max", maxdepth)}

# input files
if (ibam == "yes") {options <- paste(options, "-ibam file.a -g file.b")}
if (ibam == "no") {options <- paste(options, "-i file.a -g file.b")}

# command
command <- paste(binary, options, " > genomecoveragebed.txt")

# run
system(command)
if (file.info("genomecoveragebed.txt")$size == 0) {system("echo \"No results found\" > genomecoveragebed.txt")}