# TOOL stat-linear-modelling-designmat.R: "Linear modelling using user-defined design matrix" (Analyzes the data using linear modelling as implemented in the R package limma and users-define design-matrix. Fold changes and p-values are reported for all effects and interactions. If desired, technical replication can be blocket by settting a blocking variable.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT design.tsv: design.tsv TYPE GENERIC
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT limma.tsv: limma.tsv
# OUTPUT foldchange.tsv: foldchange.tsv 
# OUTPUT pvalues.tsv: pvalues.tsv 
# PARAMETER technical.replication: technical.replication TYPE METACOLUMN_SEL DEFAULT EMPTY (Technical replication)
# PARAMETER cont.string: contrasts TYPE STRING DEFAULT empty (List of contrasts to be compared separated by commas)
# PARAMETER p.value.adjustment.method: p.value.adjustment.method TYPE [none: none, bonferroni: bonferroni, holm: holm, hochberg: hochberg, BH: BH, BY: BY] DEFAULT BH (Multiple testing correction method)

# MK, 24.06.2013 created linear Modelling using limma
 
# Loads the libraries
library(limma)

# Loads the normalized data
file	<- c("normalized.tsv")
dat		<- read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls	<- dat[,grep("flag", names(dat))]
dat2	<- dat[,grep("chip", names(dat))]

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Loads design-matrix
design	<- read.table("design.tsv", header=T, sep="\t", row.names=1)
colnames(design) <- paste("x", colnames(design), sep="")
rownames(design) <- colnames(dat2)

# Fit linear models
if(technical.replication == "EMPTY") {
	fit				<- lmFit(dat2, design);
} else {
	techrep			<- phenodata[,grep(technical.replication, colnames(phenodata))]
	corfit			<- duplicateCorrelation(dat2, ndups=1, block=techrep)
	fit				<- lmFit(dat2, design, block=techrep, cor=corfit$consensus)
}

# Create a contrast mat and check validity
if(cont.string != "empty") {
	cont.string 	<- gsub("\\s+", "", cont.string)
	contrast.string <- paste("x", unlist(strsplit(cont.string, ",")), sep="")
	contrast.string <- gsub("-", "-x", contrast.string)
		
	for(i in 1:length(contrast.string)) {
		col.names <- unlist(strsplit(contrast.string[i], "-"))
		if(length(col.names) == 1) {
			if(!(col.names[1] %in% colnames(design))) {
				stop("CHIPSTER-NOTE: Design-matrix has no column for one of the given constrast")
			}
		} else {
			if(!(col.names[1] %in% colnames(design)) | !(col.names[2] %in% colnames(design))) {
				stop("CHIPSTER-NOTE: Design-matrix has no column for one of the given constrast")
			}
		}
	}
		
	contrasts.mat	<- makeContrasts(contrasts=contrast.string, levels=design)	
	fit2            <- contrasts.fit(fit, contrasts.mat);	
	fit				<- fit2
}

# Compute p-values
fit     			<- eBayes(fit)
colnames(fit$coef)  <- gsub("^x", "", colnames(fit$coef))
colnames(fit$coef)  <- gsub("-x", "-", colnames(fit$coef))

# Extracting data
m		<- matrix(nrow=nrow(fit$coef), ncol=ncol(fit$coef))
mm		<- matrix(nrow=nrow(fit$coef), ncol=ncol(fit$coef))
for(i in 1:ncol(fit$coef)) {
	pp			<-toptable(fit, coef=i, number=nrow(dat2), adjust.method=p.value.adjustment.method, sort.by="none")
	m[,(i)] 	<-pp$P.Value
	mm[,(i)] 	<-pp$logFC
}

# Fold change
fc					<- data.frame(mm)
fc					<- data.frame(round(fc, digits=2))
colnames(fc)		<- paste("FC.", colnames(fit$coef), sep="")
rownames(fc)		<- rownames(dat)
fc2					<- data.frame(fc)
colnames(fc2)		<- paste("chip.FC.", colnames(fit$coef), sep="")

# P-values
pvalues				<- data.frame(m)
pvalues				<- data.frame(round(pvalues, digits=6))
colnames(pvalues)	<- paste("p.adjusted.", colnames(fit$coef), sep="")
rownames(pvalues)	<- rownames(dat)
pvalues2			<- data.frame(pvalues)
colnames(pvalues2)	<- paste("chip.p.adjusted.", colnames(fit$coef), sep="")

# Create data frames
dat3				<- cbind(dat, pvalues, fc)

# Write data to disk
write.table(dat3, file="limma.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(fc2, file="foldchange.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(pvalues2, file="pvalues.tsv", sep="\t", row.names=T, col.names=T, quote=F)
