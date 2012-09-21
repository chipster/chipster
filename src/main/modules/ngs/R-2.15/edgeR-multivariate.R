# TOOL edgeR-multivariate.R: "Differential expression analysis using edgeR for multivariate experiments" (Differential expression analysis of genes from multivariate designs using the edgeR Bioconductor package. You can create the input count table and phenodata file using the tool Utilities - Define NGS experiment.)
# INPUT data.tsv TYPE GENERIC
# INPUT phenodata.tsv TYPE GENERIC
# OUTPUT OPTIONAL de-list-edger.tsv
# PARAMETER normalization: "Apply normalization" TYPE [yes, no] DEFAULT yes (Should normalization based on the trimmed mean of M-values \(TMM\) be performed to reduce the effect from sequencing biases.)
# PARAMETER main.effect1: main.effect1 TYPE METACOLUMN_SEL DEFAULT group (Main effect 1)
# PARAMETER main.effect2: main.effect2 TYPE METACOLUMN_SEL DEFAULT EMPTY (Main effect 2)
# PARAMETER main.effect3: main.effect3 TYPE METACOLUMN_SEL DEFAULT EMPTY (Main effect 3)
# PARAMETER treat.main.effect1.as.factor: treat.main.effect1.as.factor TYPE [no: no, yes: yes] DEFAULT no (Should main.effect1 be treated as a factor)
# PARAMETER treat.main.effect2.as.factor: treat.main.effect2.as.factor TYPE [no: no, yes: yes] DEFAULT no (Should main.effect2 be treated as a factor)
# PARAMETER treat.main.effect3.as.factor: treat.main.effect3.as.factor TYPE [no: no, yes: yes] DEFAULT no (Should main.effect3 be treated as a factor)
# PARAMETER interactions: interactions TYPE [main: "main effects only", all: "main effects and interactions"] DEFAULT main (What to include in the model)



# JTT 8.7.2012

# Demo settings
#normalization<-"yes"
#main.effect1<-"group"
#main.effect2<-"type"
#main.effect3<-"EMPTY"
#treat.main.effect1.as.factor<-"no"
#treat.main.effect2.as.factor<-"no"
#treat.main.effect3.as.factor<-"no"
#interactions<-"all"
#setwd("C:/Users/Jarno Tuimala/Desktop/c")


# Loads the libraries
library(edgeR)

# Loads the count data
dat <- read.table("data.tsv", header=T, sep="\t", row.names=1)

# Separates expression values and flags
annotations <- dat[,-grep("chip", names(dat))]
dat2 <- dat[,grep("chip", names(dat))]

# Reads the phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")

# Forms the DGElist object
dge<-DGEList(counts=dat2)

# Calculate normalization factors
if(normalization=="yes") {
   dge<-calcNormFactors(dge)
}

# Form a model matrix
formula<-"~"
if(main.effect1!="EMPTY" & treat.main.effect1.as.factor=="no") {
   formula<-paste(formula, main.effect1, sep="")
}
if(main.effect1!="EMPTY" & treat.main.effect1.as.factor=="yes") {
   formula<-paste(formula, "as.factor(", main.effect1, ")", sep="")
}

if(interactions=="main" & main.effect2!="EMPTY") {
   formula<-paste(formula, "+", sep="")
} 
if(interactions=="all" & main.effect2!="EMPTY") {
   formula<-paste(formula, "*", sep="")
} 
  
if(main.effect2!="EMPTY" & treat.main.effect2.as.factor=="no") {
   formula<-paste(formula, main.effect2, sep="")
}
if(main.effect2!="EMPTY" & treat.main.effect2.as.factor=="yes") {
   formula<-paste(formula, "as.factor(", main.effect2, ")", sep="")
}

if(interactions=="main" & main.effect3!="EMPTY") {
   formula<-paste(formula, "+", sep="")
} 
if(interactions=="all" & main.effect3!="EMPTY") {
   formula<-paste(formula, "*", sep="")
} 
 
if(main.effect3!="EMPTY" & treat.main.effect3.as.factor=="no") {
   formula<-paste(formula, main.effect3, sep="")
}
if(main.effect3!="EMPTY" & treat.main.effect3.as.factor=="yes") {
   formula<-paste(formula, "as.factor(", main.effect3, ")", sep="")
}

design<-with(phenodata, model.matrix(as.formula(formula)))

# Estimate disersions
dge <- estimateCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# Estimate DE genes
fit<-glmFit(dge, design)

# LRT
lrt<-glmLRT(dge, fit, coef=1)
tt<-topTags(lrt, n=nrow(dat2))
tt<-tt@.Data[[1]]
colnames(tt)<-paste(colnames(tt), colnames(design)[1], sep="-")
ttres<-tt[order(rownames(tt)),]

for(i in 2:ncol(design)) {
   lrt<-glmLRT(dge, fit, coef=i)
   tt<-topTags(lrt, n=nrow(dat2))
   tt<-tt@.Data[[1]]
   colnames(tt)<-paste(colnames(tt), colnames(design)[i], sep="-")
   tt<-tt[order(rownames(tt)),]
   ttres<-cbind(ttres, tt)
}

write.table(ttres, file="de-list-edger.tsv", sep="\t", row.names=T, col.names=T, quote=F)

