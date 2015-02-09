# TOOL acgh-survival-test.R: "Survival test for called copy number data" (Statistical test for survival and called copy number data. The testing is recommended to be performed after running the Identify common regions from called copy number data tool.)
# INPUT regions.tsv: regions.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT survival-test.tsv: survival-test.tsv 
# PARAMETER survival: Survival TYPE METACOLUMN_SEL DEFAULT survival (Phenodata column with survival data)
# PARAMETER status: Status TYPE METACOLUMN_SEL DEFAULT status (Phenodata column with patient status: alive=0, dead=1)
# PARAMETER number.of.permutations: "Number of permutations" TYPE INTEGER DEFAULT 10000 (The number of permutations. At least 10000 recommended for final calculations.)
# PARAMETER test.aberrations: "Test aberrations" TYPE [1: gains, -1: losses, 0: both] DEFAULT 0 (Whether to test only for gains or losses, or both.) 

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-23

source(file.path(chipster.common.path, 'library-Chipster.R'))

dat <- readData("regions.tsv")
phenodata <- readPhenodata("phenodata.tsv")

data.info <- dat[, annotationColumns(dat)]
calls <- as.matrix(dat[,grep('^flag\\.', colnames(dat))])

# first try parallel computing
prob <- TRUE
try({
  library(CGHtestpar)
  pvs <-  pvalstest_logrank(calls, data.info, dataclinvar=phenodata, whtime=which(colnames(phenodata) == survival), whstatus=which(colnames(phenodata) == status), lgonly=as.integer(test.aberrations), niter=number.of.permutations, ncpus=4)
  fdrs <- fdrperm(pvs)
  prob <- FALSE
}, silent=TRUE)
# if problems, fall back to sequential computing
if (prob) {
  library(CGHtest)
  pvs <-  pvalstest_logrank(calls, data.info, dataclinvar=phenodata, whtime=which(colnames(phenodata) == survival), whstatus=which(colnames(phenodata) == status), lgonly=as.integer(test.aberrations), niter=number.of.permutations)
  fdrs <- fdrperm(pvs)
}

fdrs <- cbind(fdrs, dat[, dataColumns(dat)])

writeData(fdrs, "survival-test.tsv")

# EOF
