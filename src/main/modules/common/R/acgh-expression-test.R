# TOOL acgh-expression-test.R: "Test for copy-number-induced expression changes" (Nonparametric testing for changes in expression induced by a change in DNA copy number. Before running this tool, the copy number and expression data must be first matched together using the Match copy number and expression features tool.)
# INPUT matched-cn-and-expression.tsv: matched-cn-and-expression.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT cn-induced-expression.tsv: cn-induced-expression.tsv 
# PARAMETER test.statistic: "Test statistic" TYPE [wcvm: wcvm, wmw: wmw] DEFAULT wcvm (The test statistic to use.)
# PARAMETER analysis.type: "Analysis type" TYPE [univariate: univariate, regional: regional] DEFAULT univariate (The type of the analysis.)
# PARAMETER number.of.permutations: "Number of permutations" TYPE INTEGER DEFAULT 10000 (The number of permutations used for the p-value calculation.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-23

source(file.path(chipster.common.path, "library-Chipster.R"))
source(file.path(chipster.common.path, "library-intCNGEan.R"))

# read the input files
input <- readData("matched-cn-and-expression.tsv")
phenodata <- readPhenodata("phenodata.tsv")

matched <- toIntCNGEanMatched(input)

# tune and test
tuned <- intCNGEan.tune(matched$CNdata.matched, matched$GEdata.matched, test.statistic=test.statistic)
result <- intCNGEan.test(tuned, analysis.type=analysis.type, test.statistic=test.statistic, nperm=number.of.permutations)

# for some reason, columns in the data.frame appear as lists in R-2.12.
# convert back to vectors
if (class(result[,1]) == "list")
  for (i in 1:ncol(result))
    result[,i] <- unlist(result[,i])

# format and write result table
rownames(result) <- rownames(tuned$ann)
colnames(result)[1:3] <- tolower(colnames(result)[1:3])
result$chromosome <- chromosomeToCharacter(result$chromosome)
result$probes <- rownames(tuned$datafortest)

result <- addAnnotationColumns(input, result, rows=result$probes, exclude=c("cn.start", "cn.end", "exp.probe", "exp.start", "exp.end"))

arrays <- colnames(tuned$datafortest)[(2*tuned$nosamp+1):(3*tuned$nosamp)]
colnames(tuned$datafortest) <- c(paste0("prob1.", arrays), paste0("prob2.", arrays), paste0("chip.", arrays))
result <- cbind(result, tuned$datafortest)
result <- result[order(result$adj.p),]

writeData(result, "cn-induced-expression.tsv")

# EOF
