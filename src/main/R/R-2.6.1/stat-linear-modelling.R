# ANALYSIS Statistics/"Linear modelling" (Analyzes the data using linear modelling as implemented in limma R package.
# You can have a maximum of three main effect and their interaction in the model. On top of the main effects,
# you can specify technical replication and pairing of the samples. Main effects can be fitted as such of as categorical
# variables, i.e., factors. Fold changes and p-values are reported for all effects and interactions.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT limma.tsv, limma-design.tsv
# PARAMETER column1 METACOLUMN_SEL DEFAULT group (Main effect 1)
# PARAMETER column2 METACOLUMN_SEL DEFAULT EMPTY (Main effect 2)
# PARAMETER column3 METACOLUMN_SEL DEFAULT EMPTY (Main effect 3)
# PARAMETER column4 METACOLUMN_SEL DEFAULT EMPTY (Technical replication)
# PARAMETER column5 METACOLUMN_SEL DEFAULT EMPTY (Paired samples)
# PARAMETER col1.factor [no, yes] DEFAULT no (Should column1 be treated as a factor)
# PARAMETER col2.factor [no, yes] DEFAULT no (Should column2 be treated as a factor)
# PARAMETER col3.factor [no, yes] DEFAULT no (Should column3 be treated as a factor)
# PARAMETER p.value.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER adjust.p.values [yes, no] DEFAULT yes (Should the p-values be adjusted for multiple comparisons)
# PARAMETER p.value.adjustment.method [none, bonferroni, holm, hochberg, BH, BY] DEFAULT BH (Multiple testing correction method)
# PARAMETER interactions [main, two-way, three-way] DEFAULT main (What to include in the model)
# PARAMETER significance [main, interactions] DEFAULT main (Which p-values to return)


# Linear Modelling using limma
# 
# JTT 22.10.2007

# Loads the libraries
library(limma)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Sanity checks
if(column1=="EMPTY" & column2=="EMPTY" & column3=="EMPTY") {
   stop("You need to specify at least one main effect! Please modify the setting accordingly, and rerun.")
}
if(column2=="EMPTY" | column3=="EMPTY" & interactions=="two-way" | interactions=="three-way") {
   print("Only one main effect specified with interactions! No interactions specified for the model.")
   interactions<-c("main")
}
if(interactions=="main" & significance=="interactions") {
   print("Interactions can't be tested, since only main effect specified for the model!")
   interactions<-c("main")
}

# Extracting the variables from phenodata
col1<-phenodata[,grep(column1, colnames(phenodata))]
col2<-phenodata[,grep(column2, colnames(phenodata))]
col3<-phenodata[,grep(column3, colnames(phenodata))]
col4<-phenodata[,grep(column4, colnames(phenodata))]
col5<-phenodata[,grep(column5, colnames(phenodata))]

# Converting vectors to factor, if needed
if(col1.factor=="yes") {
   col1<-factor(col1)
}
if(col2.factor=="yes") {
   col2<-factor(col2)
}
if(col3.factor=="yes") {
   col3<-factor(col3)
}

# Recoding variables
pcut<-p.value.threshold

# Specifying the models

# The basic models

# One main effect
if(column1!="EMPTY" & column2=="EMPTY" & column3=="EMPTY" & column4=="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col1)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2!="EMPTY" & column3=="EMPTY" & column4=="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col2)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2=="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col3)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}

# Two main effects
if(column1!="EMPTY" & column2!="EMPTY" & column3=="EMPTY" & column4=="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col1+col2)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}
if(column1!="EMPTY" & column2=="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col1+col3)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col2+col3)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}

# Three main effects
if(column1!="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col1+col2+col3)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}

# Two main effects and interactions
if(column1!="EMPTY" & column2!="EMPTY" & column3=="EMPTY" & column4=="EMPTY" & column5=="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col1*col2)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}
if(column1!="EMPTY" & column2=="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5=="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col1*col3)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5=="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col2*col3)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}

# Three main effects
if(column1!="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col1+col2+col3)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}

# Three main effects and interactions
if(column1!="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5=="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col1*col2*col3)
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}


# The basic models and biological replication

# Only technical replication
if(column1=="EMPTY" & column2=="EMPTY" & column3=="EMPTY" & column4!="EMPTY" & column5=="EMPTY") {
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, block=col4, cor=corfit$consensus) 
   fit<-eBayes(fit)
}

# One main effect and technical replication
if(column1!="EMPTY" & column2=="EMPTY" & column3=="EMPTY" & column4!="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col1)
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2!="EMPTY" & column3=="EMPTY" & column4!="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col2)
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2=="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col3)
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}

# Two main effects and technical replication
if(column1!="EMPTY" & column2!="EMPTY" & column3=="EMPTY" & column4!="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col1+col2)
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}
if(column1!="EMPTY" & column2=="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col1+col3)
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col2+col3)
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}

# Three main effects and technical replication
if(column1!="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col1+col2+col3)
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}

# Two main effects and interactions and technical replication
if(column1!="EMPTY" & column2!="EMPTY" & column3=="EMPTY" & column4!="EMPTY" & column5=="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col1*col2)
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}
if(column1!="EMPTY" & column2=="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5=="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col1*col3)
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5=="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col2*col3)
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}

# Three main effects and technical replication
if(column1!="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5=="EMPTY") {
   design<-model.matrix(~col1+col2+col3)
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}

# Three main effects and interactions and technical replication
if(column1!="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5=="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col1*col2*col3)
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}


# Paired models

# One main effect
if(column1!="EMPTY" & column2=="EMPTY" & column3=="EMPTY" & column4=="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col1+factor(col5))
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2!="EMPTY" & column3=="EMPTY" & column4=="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col2+factor(col5))
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2=="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col3+factor(col5))
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}

# Two main effects
if(column1!="EMPTY" & column2!="EMPTY" & column3=="EMPTY" & column4=="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col1+col2+factor(col5))
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}
if(column1!="EMPTY" & column2=="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col1+col3+factor(col5))
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col2+col3+factor(col5))
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}

# Three main effects
if(column1!="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col1+col2+col3+factor(col5))
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}

# Two main effects and interactions
if(column1!="EMPTY" & column2!="EMPTY" & column3=="EMPTY" & column4=="EMPTY" & column5!="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col1*col2+factor(col5))
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}
if(column1!="EMPTY" & column2=="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5!="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col1*col3+factor(col5))
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5!="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col2*col3+factor(col5))
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}

# Three main effects
if(column1!="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col1+col2+col3+factor(col5))
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}

# Three main effects and interactions
if(column1!="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4=="EMPTY" & column5!="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col1*col2*col3+factor(col5))
   fit<-lmFit(dat2, design)
   fit<-eBayes(fit)
}

# The basic models and biological replication

# Only technical replication
if(column1=="EMPTY" & column2=="EMPTY" & column3=="EMPTY" & column4!="EMPTY" & column5!="EMPTY") {
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, block=col4, cor=corfit$consensus) 
   fit<-eBayes(fit)
}

# One main effect and technical replication
if(column1!="EMPTY" & column2=="EMPTY" & column3=="EMPTY" & column4!="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col1+factor(col5))
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2!="EMPTY" & column3=="EMPTY" & column4!="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col2+factor(col5))
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2=="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col3+factor(col5))
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}

# Two main effects and technical replication
if(column1!="EMPTY" & column2!="EMPTY" & column3=="EMPTY" & column4!="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col1+col2+factor(col5))
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}
if(column1!="EMPTY" & column2=="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col1+col3+factor(col5))
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col2+col3+factor(col5))
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}

# Three main effects and technical replication
if(column1!="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col1+col2+col3+factor(col5))
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}

# Two main effects and interactions and technical replication
if(column1!="EMPTY" & column2!="EMPTY" & column3=="EMPTY" & column4!="EMPTY" & column5!="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col1*col2+factor(col5))
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}
if(column1!="EMPTY" & column2=="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5!="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col1*col3+factor(col5))
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}
if(column1=="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5!="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col2*col3+factor(col5))
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}

# Three main effects and technical replication
if(column1!="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5!="EMPTY") {
   design<-model.matrix(~col1+col2+col3+factor(col5))
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}

# Three main effects and interactions and technical replication
if(column1!="EMPTY" & column2!="EMPTY" & column3!="EMPTY" & column4!="EMPTY" & column5!="EMPTY" & interactions=="two-way") {
   design<-model.matrix(~col1*col2*col3+factor(col5))
   corfit<-duplicateCorrelation(dat2, ndups=1, block=col4)
   fit<-lmFit(dat2, design, block=col4, cor=corfit$consensus)
   fit<-eBayes(fit)
}


# Extracting the significant results
# m<-matrix(nrow=nrow(dat2), ncol=ncol(design)-1)
# mm<-matrix(nrow=nrow(dat2), ncol=ncol(design)-1)
# for(i in 2:ncol(design)) {
#   pp<-toptable(fit, coef=i, number=nrow(dat2), adjust.method=p.value.adjustment.method)
#   ppind<-as.numeric(rownames(pp))
#   ppind<-order(ppind)
#   if(adjust.p.values=="yes") {
#      pp2<-pp$adj.P.Val[ppind]
#   }
#   if(adjust.p.values=="no") {
#      pp2<-pp$P.Value[ppind]
#   }
#   m[,(i-1)]<-pp2
#   pp3<-pp$logFC[ppind]
#   mm[,(i-1)]<-pp3
#}
#
# Writing the data to disk
#m<-data.frame(m)
#mm<-data.frame(mm)
#colnames(m)<-paste("p.adjusted.", colnames(design)[-1], sep="")
#colnames(mm)<-paste("FC.", colnames(design)[-1], sep="")
#dat3<-data.frame(dat2, round(m, digits=4), round(mm, digits=2))
#write.table(dat3, file="limma.tsv", sep="\t", row.names=T, col.names=T, quote=F)
#write.table(design, file="limma-design.tsv", sep="\t", row.names=F, col.names=T, quote=F)


## New method ##
m<-matrix(nrow=nrow(dat2), ncol=ncol(design))
mm<-matrix(nrow=nrow(dat2), ncol=ncol(design))
for(i in 1:ncol(design)) {
   pp<-toptable(fit, coef=i, number=nrow(dat2), adjust.method=p.value.adjustment.method, sort.by="p")
   rownames(pp)<-rownames(dat2[as.numeric(rownames(pp)),])
   if(adjust.p.values=="yes") {
      pp2<-pp$adj.P.Val
   }
   if(adjust.p.values=="no") {
      pp2<-pp$P.Value
   }
   pp3<-pp$logFC
   pp2<-pp2[order(rownames(pp))]
   pp3<-pp3[order(rownames(pp))]
   m[,(i)]<-pp2
   mm[,(i)]<-pp3
}

# Writing the data to disk
m<-data.frame(m)
mm<-data.frame(mm)
colnames(m)<-paste("p.adjusted.", colnames(design), sep="")
colnames(mm)<-paste("FC.", colnames(design), sep="")
dat3<-data.frame(dat, round(m, digits=4), round(mm, digits=2))
write.table(dat3, file="limma.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(design, file="limma-design.tsv", sep="\t", row.names=F, col.names=T, quote=F)
