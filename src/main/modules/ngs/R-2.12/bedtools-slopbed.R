# TOOL bedtools-slopbed.R: "Slop BED" (Increases the size of each feature in a feature file be a user-defined number of bases, but restricts the resizing to the size of the chromosome (i.e. no start < 0 and no end > chromosome size\). This tool is based on the BEDTools package.)
# INPUT bed.file: "BED/GFF/VCF file" TYPE GENERIC
# INPUT gen.file: "Genome file" TYPE GENERIC
# OUTPUT slopbed.bed
# PARAMETER OPTIONAL l: "The number of base pairs to subtract from the start coordinate" TYPE INTEGER DEFAULT 0 ()
# PARAMETER OPTIONAL r: "The number of base pairs to add to the end coordinate" TYPE INTEGER DEFAULT 0 ()
# PARAMETER OPTIONAL s: "Define additions based on strand" TYPE [yes, no] DEFAULT no (If this option is selected and for example 500 bases are subtracted from a start coordinate of a negative-stranded feature,  500 bp will be added to the end coordinate.)
# PARAMETER OPTIONAL pct: "Define additions as a fraction of the feature's length." TYPE [yes, no] DEFAULT no (Define additions as a fraction of the feature's length. E.g. if used on a 1000bp feature, 0.50 will add 500 bp.)
# PARAMETER OPTIONAL lf: "The fraction to subtract from the start coordinate" TYPE DECIMAL DEFAULT 0 ()
# PARAMETER OPTIONAL rf: "The fraction to add to the end coordinate" TYPE DECIMAL DEFAULT 0 ()

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "slopBed"))

# options
options <- paste("")
if (pct == "no") {options <- paste(options, "-l", l, "-r", r)}
if (pct == "yes") {options <- paste(options, "-pct", "-l", lf, "-r", rf)}
if (s == "yes") {options <- paste(options, "-s")}

# input files
options <- paste(options, "-i bed.file", "-g gen.file")



# command
command <- paste(binary, options, "> slopbed.bed")

# run
system(command)
if (file.info("slopbed.bed")$size == 0) {system("echo \"No results found\" > slopbed.bed")}