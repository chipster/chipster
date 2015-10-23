# TOOL edgeR-multivariate.R: "Differential expression using edgeR for multivariate experiments" (Differential expression analysis for multifactor experiments using the generalized linear models based statistical methods of the edgeR Bioconductor package. You can create the input count table and phenodata file using the tool "\Utilities - Define NGS experiment\".)
# INPUT data.tsv TYPE GENERIC
# INPUT META phenodata.tsv TYPE GENERIC
# OUTPUT OPTIONAL edger-glm.tsv
# OUTPUT OPTIONAL dispersion-edger-glm.pdf
# PARAMETER main.effect1: "Main effect 1" TYPE METACOLUMN_SEL DEFAULT group (Main effect 1)
# PARAMETER OPTIONAL main.effect2: "Main effect 2" TYPE METACOLUMN_SEL DEFAULT EMPTY (Main effect 2)
# PARAMETER OPTIONAL main.effect3: "Main effect 3" TYPE METACOLUMN_SEL DEFAULT EMPTY (Main effect 3)
# PARAMETER OPTIONAL treat.main.effect1.as.factor: "Treat main effect 1 as factor" TYPE [no: no, yes: yes] DEFAULT yes (Should main.effect1 be treated as a factor)
# PARAMETER OPTIONAL treat.main.effect2.as.factor: "Treat main effect 2 as factor" TYPE [no: no, yes: yes] DEFAULT yes (Should main.effect2 be treated as a factor)
# PARAMETER OPTIONAL treat.main.effect3.as.factor: "Treat main effect 3 as factor" TYPE [no: no, yes: yes] DEFAULT yes (Should main.effect3 be treated as a factor)
# PARAMETER OPTIONAL interactions: "Include interactions in the model" TYPE [main: "no", all: "yes", nested: "nested"] DEFAULT main (Should interactions be included in the model in addition to the main effects. "Yes" = include all interaction terms, "no" = include only main effects, "nested" = include main effect 1 and its interactions with effects 2 and 3. Use the last option when you have comparisons both between and within subjects. Please note that when using the "nested" option, you need to place the effects in particular order -see manual for more information! )
# PARAMETER OPTIONAL normalization: "Apply TMM normalization" TYPE [yes, no] DEFAULT yes (Should normalization based on the trimmed mean of M-values \(TMM\) be performed to reduce the RNA composition effect.)
# PARAMETER OPTIONAL filter: "Analyze only genes which have counts in at least this many samples" TYPE INTEGER FROM 0 TO 1000 DEFAULT 1 (Analyze only genes which have at least 5 counts in at least this many samples)
# PARAMETER OPTIONAL w: "Plot width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted image)
# PARAMETER OPTIONAL h: "Plot height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted image)


# JTT 8.7.2012
# EK 28.4.2013 rounding added, main effect treated as factor by default
# EK 2.5.2013 updated to BioC2.11
# EK 4.5.2013 added dispersion plot and filtering
# EK 19.11.2013 updated to edgeR 3.4.0 and added the counts to output
# ML & SS 3.3.2015 added the "nested" option to the interactions


# Loads the libraries
library(edgeR)

# Loads the count data
dat <- read.table("data.tsv", header=T, sep="\t", row.names=1)

# Extracts expression value columns
annotations <- dat[,-grep("chip", names(dat))]
dat2 <- dat[,grep("chip", names(dat))]

# Reads the phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")

# Forms the DGElist object
dge<-DGEList(counts=dat2)

# filter out genes which have less than 5 counts in user-defined number of samples 
if (filter > 0) {
	keep <- rowSums(dge$counts>5) >= filter
	dge <- dge[keep,]
	dge$lib.size <- colSums(dge$counts)
}

# Calculate normalization factors
if(normalization=="yes") {
	dge<-calcNormFactors(dge)
}
# For later use: filter out genes which have less than 1 cpm in user-defined number of samples 
# if (filter > 0) {
# 	keep <- rowSums(cpm(dge)>1) >= filter
# 	dge <- dge[keep,]
# 	dge$lib.size <- colSums(dge$counts)
# }

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

# If there are interactions between and within the subjects (see edgeR userguide chapter 3.5)
if(main.effect3!="EMPTY" & interactions=="nested") {
	formula<-paste("~as.factor(", main.effect1,")+as.factor(", main.effect1, "):as.factor(", main.effect2,")+as.factor(", main.effect1, "):as.factor(", main.effect3, ")", sep="")

	# modify the phenodata column for "subject", create a tmp table for this:
	tmp <- phenodata[c(main.effect1, main.effect2, main.effect3)]
	tmp[[main.effect1]]<-factor(tmp[[main.effect1]]) # group main effect 1, (disease group,) "between" comparison
	tmp[[main.effect2]]<-factor(tmp[[main.effect2]]) # trt main effect 2, (treatment,) the "within" comparison
	tmp[[main.effect3]]<-factor(tmp[[main.effect3]])  # subj main effect 3, (patient,) pairing
	
	splitsubj<-split(tmp[[main.effect3]],tmp[[main.effect1]])
	splitind<-lapply(splitsubj,FUN=function(ssubj) {as.numeric(factor(ssubj))} )
	## in the anonymous function, call to factor will drop unused levels, 
	## and as.numeric returns only the integer codes (so the result is not a factor) to problems in unsplit
	tmp[[main.effect3]]<-unsplit(splitind,tmp[[main.effect1]])
	
	phenodata <- tmp

}

design<-with(phenodata, model.matrix(as.formula(formula)))
# remove zero columns:
design <- design[, colSums(abs(design),na.rm = TRUE) != 0]

# Estimate dispersions
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# Dispersion plot
pdf(file="dispersion-edger-glm.pdf", width=w/72, height=h/72)
plotBCV(dge, main="Biological coefficient of variation")
dev.off()

# Estimate DE genes
fit<-glmFit(dge, design)

# LRT
lrt<-glmLRT(fit, coef=1)
tt<-topTags(lrt, n=nrow(dat2))
tt<-tt@.Data[[1]]
colnames(tt)<-paste(colnames(tt), colnames(design)[1], sep="-")
ttres<-tt[order(rownames(tt)),]

for(i in 2:ncol(design)) {
	lrt<-glmLRT(fit, coef=i)
	tt<-topTags(lrt, n=nrow(dat2))
	tt<-tt@.Data[[1]]
	colnames(tt)<-paste(colnames(tt), colnames(design)[i], sep="-")
	tt<-tt[order(rownames(tt)),]
	ttres<-cbind(ttres, tt)
}

# Rounding, etc.
ttres2<-round(ttres,6)

# add count columns to the result table
ttres2<-cbind(ttres2, getCounts(dge))


#ttres3<-merge(dat, ttres2, by.x=0, by.y=0)

write.table(ttres2, file="edger-glm.tsv", sep="\t", row.names=T, col.names=T, quote=F)

