# TOOL vcftools-reports.R: "Calculate statistics on VCF file" (Given a VCF file, this tool calculates several statistics such as allele frequency and count, LD measure, and p-value for HWE. It can also produce the number and density of SNPs in bins of defined size. This tool is based on the VCFtools package.)
# INPUT input.vcf: "VCF file" TYPE GENERIC 
# OUTPUT OPTIONAL vcftools.log
# OUTPUT OPTIONAL vcftools.frq.tsv
# OUTPUT OPTIONAL vcftools.frq.count.tsv
# OUTPUT OPTIONAL vcftools.hwe.tsv
# OUTPUT OPTIONAL vcftools.geno.ld.tsv
# OUTPUT OPTIONAL vcftools.snpden.tsv
# PARAMETER OPTIONAL statistics.freq: "Report per-site frequency information" TYPE [yes, no] DEFAULT no (Reports per-site allele frequency and count information.)
# PARAMETER OPTIONAL statistics.pvalue: "Report p-value" TYPE [yes, no] DEFAULT no (Reports a p-value for each site from a Hardy-Weinberg Equilibrium test.)
# PARAMETER OPTIONAL statistics.ld: "Report LD statistics" TYPE [yes, no] DEFAULT no (Report Linkage Disequilibrium (LD\) statistics.)
# PARAMETER OPTIONAL statistics.snpdensity: "Report SNP density" TYPE INTEGER DEFAULT 0 (Calculates the number and density of SNPs in bins of size given by the user. If the size is 0, this option is ignored.)

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

# statistics options

if (statistics.freq == "yes"){
	command <- paste(vcftools.binary, vcftools.io.options, "--freq", "2>> vcftools.log")
	system(command)
	command <- paste(vcftools.binary, vcftools.io.options, "--counts", "2>> vcftools.log")
	system(command)
}
if (statistics.pvalue == "yes"){
	command <- paste(vcftools.binary, vcftools.io.options, "--hardy", "2>> vcftools.log")
	system(command)
}
if (statistics.ld == "yes"){
	command <- paste(vcftools.binary, vcftools.io.options, "--geno-r2", "2>> vcftools.log")
	system(command)
}
if (statistics.snpdensity != "0"){
	command <- paste(vcftools.binary, vcftools.io.options, "--SNPdensity", statistics.snpdensity, "2>> vcftools.log")
	system(command)
}

# Fix last column on .frq and .frq.count to be Chipster tsv compatible
awk.command <- paste("| awk \'{printf $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5; for(i=6;i<=NF;i++){printf \",%s\", $i} printf \"\\n\"}\'")
system(paste("cat vcftools.frq", awk.command, "> vcftools.frq.tsv"))
system(paste("cat vcftools.frq.count", awk.command, "> vcftools.frq.count.tsv"))

# rename result files
system("mv vcftools.hwe vcftools.hwe.tsv")
system("mv vcftools.geno.ld vcftools.geno.ld.tsv")
system("mv vcftools.snpden vcftools.snpden.tsv")
