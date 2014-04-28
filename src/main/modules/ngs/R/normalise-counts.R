# TOOL normalise-counts.R: "Normalise count table" (Converts count values to counts-per-million values. Conversion is done using functions from edgeR package and can include TMM, RLE, upperquantile normalisation.)
# INPUT counts.tsv: "RNA count data" TYPE GENE_EXPRS 
# OUTPUT OPTIONAL ngs-data-table.tsv: ngs-data-table.tsv
# PARAMETER method.value: "Calculate normalization factors" TYPE [TMM: TMM, RLE: RLE, upperquantile: upperquantile, none: none] DEFAULT TMM (Method to be used for computing normalization factors. None corresponds to lib.size normalisation) 
# PARAMETER log.value: "Log transform" TYPE [TRUE: yes, FALSE: no] DEFAULT yes (Is data log2-transformed.) 

# 25.04.2014 MK, created

# Loads the normalized data and phenodata files
data <- read.table(file="counts.tsv", header=T, sep="\t", row.names=1)
count.data <- data[,grep("chip", names(data))]

library(edgeR)
cds <- DGEList(count.data)
cds <- calcNormFactors(cds, method=method.value)

# Normalization
if(method.value=="none") {
  data2 <- cpm(cds, normalized.lib.sizes=FALSE, log=log.value)
} else {
  data2 <- cpm(cds, normalized.lib.sizes=TRUE, log=log.value)  
}

# Write out results
data[,grep("chip", names(data))] <- data2
write.table(data, "ngs-data-table.tsv", col.names=T, row.names=T, sep="\t", quote=F)