# TOOL bcftools-intersect-vcf.R: "Intersect VCF" (Makes an intersection of two VCF files. This tool is based on the BCFtools package.)
# INPUT a.vcf: "VCF file A" TYPE GENERIC
# INPUT b.vcf: "VCF file B" TYPE GENERIC
# OUTPUT OPTIONAL 0000.vcf
# OUTPUT OPTIONAL 0001.vcf
# OUTPUT OPTIONAL 0002.vcf
# OUTPUT OPTIONAL 0003.vcf

# AMS 02.10.2015

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("a.vcf")
unzipIfGZipFile("b.vcf")

# binaries
bcftools.binary <- c(file.path(chipster.tools.path, "bcftools", "bcftools"))
bgzip.binary <- c(file.path(chipster.tools.path, "tabix", "bgzip"))

# Bgzip and index input files
system(paste(bgzip.binary, "a.vcf"))
system(paste(bgzip.binary, "b.vcf"))
system(paste(bcftools.binary, "index", "a.vcf.gz"))
system(paste(bcftools.binary, "index", "b.vcf.gz"))

# Intersect VCFs
system(paste(bcftools.binary, " isec -O v -p tmp a.vcf.gz b.vcf.gz"))

# rename result files
system("mv tmp/0000.vcf 0000.vcf")
system("mv tmp/0001.vcf 0001.vcf")
system("mv tmp/0002.vcf 0002.vcf")
system("mv tmp/0003.vcf 0003.vcf")

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))

# Make a matrix of output names
outputnames <- matrix(NA, nrow=4, ncol=2)
outputnames[1,] <- c("0000.vcf", paste("variantsOnlyInA.vcf"))
outputnames[2,] <- c("0001.vcf", paste("variantsOnlyInB.vcf"))
outputnames[3,] <- c("0002.vcf", paste("variantsInA_matchingB.vcf"))
outputnames[4,] <- c("0003.vcf", paste("variantsInB_matchingA.vcf"))

# Write output definitions file
write_output_definitions(outputnames)
