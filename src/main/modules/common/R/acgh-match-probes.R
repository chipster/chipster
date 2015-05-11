# TOOL acgh-match-probes.R: "Match copy number and expression features" (Matches the data points of a copy number data set with data points of an expression data set using their chromosomal locations. Running this tool is a prerequisite for testing copy-number-induced effects on expression. See manual page for more details.)
# INPUT aberrations.tsv: "copy number data" TYPE GENE_EXPRS 
# INPUT normalized.tsv: "expression data" TYPE GENE_EXPRS 
# INPUT META phenodata_cgh.tsv: "copy number phenodata" TYPE GENERIC 
# INPUT META phenodata_exp.tsv: "expression phenodata" TYPE GENERIC 
# OUTPUT matched-cn-and-expression.tsv: matched-cn-and-expression.tsv 
# OUTPUT matched-cn-and-expression-heatmap.pdf: matched-cn-and-expression-heatmap.pdf 
# OUTPUT META phenodata-matched.tsv: phenodata-matched.tsv 
# PARAMETER samples_cgh: "Copy number sample identifiers" TYPE METACOLUMN_SEL DEFAULT Sample (The phenodata column for data set 1 used to link the two data sets together.)
# PARAMETER samples_exp: "Expression sample identifiers" TYPE METACOLUMN_SEL DEFAULT Sample (The phenodata column for data set 2 used to link the two data sets together.)
# PARAMETER method: Method TYPE [distance: distance, overlap: overlap, overlapplus: overlapplus] DEFAULT distance (The method for linking copy number and expression data points together.)

# match-cn-and-expression-probes.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-25

source(file.path(chipster.common.path, "library-Chipster.R"))
source(file.path(chipster.common.path, "library-CGHcall.R"))
source(file.path(chipster.common.path, "library-ExpressionSet.R"))
source(file.path(chipster.common.path, "library-intCNGEan.R"))

input_cgh <- readData("aberrations.tsv")
phenodata_cgh <- readPhenodata("phenodata_cgh.tsv")

input_exp <- readData("normalized.tsv")
phenodata_exp <- readPhenodata("phenodata_exp.tsv")

# check for unambiguity of sample identifiers
if (nrow(phenodata_cgh)!=length(unique(phenodata_cgh[,samples_cgh])))
  stop("CHIPSTER-NOTE: Unambigous copy number sample identifiers: ", paste(phenodata_cgh[duplicated(phenodata_cgh[,samples_cgh]),samples_cgh], collapse=", ")) 
if (nrow(phenodata_exp)!=length(unique(phenodata_exp[,samples_exp])))
  stop("CHIPSTER-NOTE: Unambigous expression sample identifiers: ", paste(phenodata_exp[duplicated(phenodata_exp[,samples_exp]),samples_exp], collapse=", ")) 

common.samples <- intersect(phenodata_cgh[,samples_cgh], phenodata_exp[,samples_exp])
rownames(phenodata_cgh) <- phenodata_cgh[,samples_cgh]
rownames(phenodata_exp) <- phenodata_exp[,samples_exp]
phenodata_cgh$n <- 1:nrow(phenodata_cgh)
phenodata_exp$n <- 1:nrow(phenodata_exp)
index_cgh <- phenodata_cgh[common.samples, "n"]
index_exp <- phenodata_exp[common.samples, "n"]

cgh <- toCgh(input_cgh, level="calls", maxStates=3)
cgh <- cgh[, index_cgh]
sampleNames(cgh) <- phenodata_cgh[common.samples, samples_cgh]

# build an ExpressionSet object from the expression data
exp <- toExpressionSet(input_exp, requirePositions=TRUE)
exp <- exp[, index_exp]
sampleNames(exp) <- phenodata_exp[common.samples, samples_exp]

# match probes
matched <- intCNGEan.match(cgh, exp, CNbpend="yes", GEbpend="yes", method=method)

# add additional phenodata columns from phenodata_cgh to phenodata_exp
phenodata_cgh <- phenodata_cgh[index_cgh,]
phenodata_exp <- phenodata_exp[index_exp,]
phenodata_exp$description_cgh <- phenodata_cgh$description
for (col in setdiff(colnames(phenodata_cgh), colnames(phenodata_exp)))
  phenodata_exp[,col] <- phenodata_cgh[,col]
phenodata_exp$sample <- sprintf("microarray%.3i", 1:row(phenodata_exp))
phenodata_exp$n <- NULL

# plot heatmaps
pdf(file="matched-cn-and-expression-heatmap.pdf", paper="a4r", width=0, height=0)
intCNGEan.heatmaps(matched$CNdata.matched, matched$GEdata.matched, location="median")
dev.off()

output <- fromIntCNGEanMatched(matched)

writeData(output, "matched-cn-and-expression.tsv")
writePhenodata(phenodata_exp, "phenodata-matched.tsv")

# EOF
