# TOOL vcftools.R: "Analyze VCF file" (Filters and analyses VCF files. This tool is based on the VCFtools package.)
# INPUT input.vcf: "VCF file" TYPE GENERIC 
# OUTPUT kept.vcf
# OUTPUT removed.vcf
# PARAMETER OPTIONAL filter.keepindels: "Keep only indels" TYPE [yes, no] DEFAULT no (Keep only indels.)
# PARAMETER OPTIONAL filter.removeindels: "Remove indels" TYPE [yes, no] DEFAULT no (Remove indels.)
# PARAMETER OPTIONAL filter.numalleles: "Filter by allele number" TYPE [yes, no] DEFAULT no (Include only sites with a number of alleles within the specified range.)
# PARAMETER OPTIONAL filter.minalleles: "Minumun number of alleles" TYPE INTEGER (Minumun number of alleles)
# PARAMETER OPTIONAL filter.maxalleles: "Maximum number of alleles" TYPE INTEGER (Maximum number of alleles)
# PARAMETER OPTIONAL statistics.freq: "Report per-site frequency information" TYPE [yes, no] DEFAULT no (Reports per-site frequency information.)
# PARAMETER OPTIONAL statistics.pvalue: "Report p-value" TYPE [yes, no] DEFAULT no (Reports a p-value for each site from a Hardy-Weinberg Equilibrium test.)
# PARAMETER OPTIONAL statistics.ld: "Report LD statistics" TYPE [yes, no] DEFAULT no (Report Linkage Disequilibrium (LD\) statistics.)
# PARAMETER OPTIONAL statistics.snpdensity: "Report SNP statistics" TYPE [yes, no] DEFAULT no (Calculates the number and density of SNPs in bins of size.)

# AMS 14.9.2012

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("input.vcf")

# binaries
vcftools.binary <- c(file.path(chipster.tools.path, "vcftools", "bin", "vcftools"))


# filtering otions
if (filter.keepindels == "yes"){
	vcftools.options <- paste(vcftools.options, "--keep-only-indels")
}
if (filter.removeindels == "yes"){
	vcftools.options <- paste(vcftools.options, "--remove-indels")
}
if (filter.numalleles == "yes"){
	vcftools.options <- paste(vcftools.options, "--min-alleles", filter.minalleles, "--max-alleles", filter.maxalleles)
}

# statistics options
if (statistics.freq == "yes"){
	vcftools.options <- paste(vcftools.options, "--freq", "--count")
}
if (statistics.pvalue == "yes"){
	vcftools.options <- paste(vcftools.options, "--hardy")
}
if (statistics.ld == "yes"){
	vcftools.options <- paste(vcftools.options, "--geno-r2")
}
if (statistics.snpdensity == "yes"){
	vcftools.options <- paste(vcftools.options, "--SNPdensity")
}


# commands
command1 <- paste(vcftools.binary, vcftools.options, "--filtered-sites")

# run
#stop(paste('CHIPSTER-NOTE: ', command1))
system(command1)

