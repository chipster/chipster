# TOOL vcftools.R: "Filter variants" (Filters variants in VCF files. This tool is based on the VCFtools package.)
# INPUT input.vcf: "VCF file" TYPE GENERIC 
# OUTPUT OPTIONAL vcftools.log
# OUTPUT OPTIONAL filtered.vcf
# OUTPUT OPTIONAL removed.sites.tsv
# PARAMETER OPTIONAL filter.keepindels: "Keep only indels" TYPE [yes, no] DEFAULT no (Keep only indels.)
# PARAMETER OPTIONAL filter.removeindels: "Remove indels" TYPE [yes, no] DEFAULT no (Remove indels.)
# PARAMETER OPTIONAL filter.minalleles: "Minimum number of alleles" TYPE INTEGER DEFAULT 0 (Minumun number of alleles. 0 value means option is ignored.)
# PARAMETER OPTIONAL filter.maxalleles: "Maximum number of alleles" TYPE INTEGER DEFAULT 0 (Maximum number of alleles. 0 value means option is ignored.)
# PARAMETER OPTIONAL filter.minquality: "Minimum quality" TYPE DECIMAL DEFAULT 0 (Include only sites with quality above this threshold. 0 value means option is ignored.)
# PARAMETER OPTIONAL filter.flag: "Filter by flags" TYPE [yes, no] DEFAULT no (Removes all sites with a FILTER flag other than PASS.)
# PARAMETER OPTIONAL output.filtered: "List removed sites" TYPE [yes, no] DEFAULT no (Creates a file listing sites that have been removed after filtering. Default is to list kept files only.)

# AMS 14.9.2012
# AMS 02.12.2012 Added option for minimum quality
# AMS 23.09.2015 Changed to be VCFtools 0.1.13 compliant (one output option only)
# AMS 24.09.2015 Split tool to two tools

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("input.vcf")

# binaries
vcftools.binary <- c(file.path(chipster.tools.path, "vcftools", "bin", "vcftools"))

# Input and output options
vcftools.io.options <- paste("--vcf input.vcf", "--out vcftools")

# filtering options
vcftools.filter.options <- paste("")
if (filter.keepindels == "yes"){
	vcftools.filter.options <- paste(vcftools.filter.options, "--keep-only-indels")
}
if (filter.removeindels == "yes"){
	vcftools.filter.options <- paste(vcftools.filter.options, "--remove-indels")
}
if (filter.minalleles != "0"){
	vcftools.filter.options <- paste(vcftools.filter.options, "--min-alleles", filter.minalleles)
}
if (filter.maxalleles != "0"){
	vcftools.filter.options <- paste(vcftools.filter.options, "--max-alleles", filter.maxalleles)
}
if (filter.minquality > 0){
	vcftools.filter.options <- paste(vcftools.filter.options, "--minQ", filter.minquality)
}
if (filter.flag == "yes"){
	vcftools.filter.options <- paste(vcftools.filter.options, "--remove-filtered-all")
}

# command
# command <- paste(vcftools.binary, vcftools.io.options, vcftools.filter.options, "--recode", "--recode-INFO-all", "2> vcftools.log")
command <- paste(vcftools.binary, vcftools.io.options, vcftools.filter.options, "--recode", "--recode-INFO DP", "2> vcftools.log")
system(command)

if (output.filtered == "yes"){
	command <- paste(vcftools.binary, vcftools.io.options, vcftools.filter.options, "--removed-sites", "2>> vcftools.log")
	system(command)
}

# rename result files
system("mv vcftools.recode.vcf filtered.vcf")
system("mv vcftools.removed.sites removed.sites.tsv")
