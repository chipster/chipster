# TOOL vcftools.R: "Filter and analyse variants" (Filters and analyses variants in VCF files. This tool is based on the VCFtools package.)
# INPUT input.vcf: "VCF file" TYPE GENERIC 
# OUTPUT vcftools.log
# OUTPUT OPTIONAL vcftools.filtered.vcf
# OUTPUT OPTIONAL vcftools.frq.tsv
# OUTPUT OPTIONAL vcftools.frq.count.tsv
# OUTPUT OPTIONAL vcftools.hwe.tsv
# OUTPUT OPTIONAL vcftools.geno.ld.tsv
# OUTPUT OPTIONAL vcftools.snpden.tsv
# OUTPUT OPTIONAL vcftools.kept.sites.tsv
# OUTPUT OPTIONAL vcftools.removed.sites.tsv
# PARAMETER OPTIONAL filter.keepindels: "Keep only indels" TYPE [yes, no] DEFAULT no (Keep only indels.)
# PARAMETER OPTIONAL filter.removeindels: "Remove indels" TYPE [yes, no] DEFAULT no (Remove indels.)
# PARAMETER OPTIONAL filter.minalleles: "Minimum number of alleles" TYPE INTEGER DEFAULT 0 (Minumun number of alleles. 0 value means option is ignored.)
# PARAMETER OPTIONAL filter.maxalleles: "Maximum number of alleles" TYPE INTEGER DEFAULT 0 (Maximum number of alleles. 0 value means option is ignored.)
# PARAMETER OPTIONAL statistics.freq: "Report per-site frequency information" TYPE [yes, no] DEFAULT no (Reports per-site frequency information.)
# PARAMETER OPTIONAL statistics.pvalue: "Report p-value" TYPE [yes, no] DEFAULT no (Reports a p-value for each site from a Hardy-Weinberg Equilibrium test.)
# PARAMETER OPTIONAL statistics.ld: "Report LD statistics" TYPE [yes, no] DEFAULT no (Report Linkage Disequilibrium (LD\) statistics.)
# PARAMETER OPTIONAL statistics.snpdensity: "Report SNP density" TYPE INTEGER DEFAULT 0 (Calculates the number and density of SNPs in bins of size. 0 value means option is ignored.)
# PARAMETER OPTIONAL output.filtered: "List removed files" TYPE [yes, no] DEFAULT no (Creates two files listing sites that have been kept or removed after filtering. Default is to list kept files only.)

# AMS 14.9.2012

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("input.vcf")

# binaries
vcftools.binary <- c(file.path(chipster.tools.path, "vcftools", "bin", "vcftools"))


# filtering otions
vcftools.options <- paste("")
if (filter.keepindels == "yes"){
	vcftools.options <- paste(vcftools.options, "--keep-only-indels")
}
if (filter.removeindels == "yes"){
	vcftools.options <- paste(vcftools.options, "--remove-indels")
}
if (filter.minalleles != "0"){
	vcftools.options <- paste(vcftools.options, "--min-alleles", filter.minalleles)
}
if (filter.maxalleles != "0"){
	vcftools.options <- paste(vcftools.options, "--max-alleles", filter.maxalleles)
}


# statistics options
if (statistics.freq == "yes"){
	vcftools.options <- paste(vcftools.options, "--freq", "--counts")
}
if (statistics.pvalue == "yes"){
	vcftools.options <- paste(vcftools.options, "--hardy")
}
if (statistics.ld == "yes"){
	vcftools.options <- paste(vcftools.options, "--geno-r2")
}
if (statistics.snpdensity != "0"){
	vcftools.options <- paste(vcftools.options, "--SNPdensity", statistics.snpdensity)
}
if (output.filtered == "yes"){
	vcftools.options <- paste(vcftools.options, "--filtered-sites")
}

# commands
command1 <- paste(vcftools.binary, "--vcf input.vcf", "--out vcftools", "--recode", vcftools.options)

# run
#stop(paste('CHIPSTER-NOTE: ', command1))
system(command1)

# rename result files
system("mv vcftools.recode.vcf vcftools.filtered.vcf")
system("mv vcftools.frq vcftools.frq.tsv")
system("mv vcftools.frq.count vcftools.frq.count.tsv")
system("mv vcftools.hwe vcftools.hwe.tsv")
system("mv vcftools.geno.ld vcftools.geno.ld.tsv")
system("mv vcftools.snpden vcftools.snpden.tsv")
system("mv vcftools.kept.sites vcftools.kept.sites.tsv")
system("mv vcftools.removed.sites vcftools.removed.sites.tsv")