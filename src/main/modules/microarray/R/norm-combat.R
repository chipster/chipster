# TOOL norm-combat.R: "ComBat - batch normalisation" (Batch normalisation using Phenotype file. You may use this method after the ordinary normalisation to remove batch effects from the data.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT batch_corrected.tsv: batch_corrected.tsv
# PARAMETER batch1: batch1 TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing the main batch to correct)
# PARAMETER batch2: batch2 TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing addittional batches to correct.)
# PARAMETER batch3: batch3 TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing addittional batches to correct.)

# ComBat batch correction analysis
# MK 25.06.2013

# Loading libraries
library(sva)

# Loads the normalized data
file <- "normalized.tsv";
dat <- read.table(file, header=T, sep="\t", row.names=1)

# Loads Batch information
file <- "phenodata.tsv";
phenodata <- read.table(file, header=T, sep="\t")

# Merge batches to matrix
if(batch1 == "EMPTY" & batch2 == "EMPTY" & batch3 == "EMPTY") { stop("CHIPSTER-NOTE: no batches defined"); }

batch <- NULL;
if(batch1 != "EMPTY") { 
	batch <- cbind(batch, phenodata[,pmatch(batch1,colnames(phenodata))]);
}
if(batch2 != "EMPTY") { 
	batch <- cbind(batch, phenodata[,pmatch(batch2,colnames(phenodata))]);
}
if(batch3 != "EMPTY") { 
	batch <- cbind(batch, phenodata[,pmatch(batch3,colnames(phenodata))]);
}

dat2  <- as.matrix(dat[,grep("chip", names(dat))])
combat_edata = ComBat(dat=dat2, batch=batch, mod=NULL)

#pheno = pData(bladderEset)
#batch = pheno$batch
#mod = model.matrix(~as.factor(cancer), data=pheno)

## Correct for batch using ComBat

cols <- grep("chip", names(dat));
dat[, cols] <- as.matrix(combat_edata);
write.table(data.frame(dat), file="batch_corrected.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

