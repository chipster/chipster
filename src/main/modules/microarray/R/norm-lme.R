# TOOL norm-lme.R: "Random effects" (Removes the possible effect of random effects from the normalized expression values of every gene. This can be used for removing batch effects due to day, technician, etc. from the data. You should have a column in your phenodata that defines the random effects groups, currently only two groups are allowed.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT normalized-lme.tsv: normalized-lme.tsv 
# PARAMETER column.groups: "Column groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column containing group effects)
# PARAMETER column.random: "Column random" TYPE METACOLUMN_SEL DEFAULT random (Phenodata column containing random effects groups)
# PARAMETER error.handle: "Error handling method" TYPE [remove: remove, na: NA, add.noise: add.noise] DEFAULT remove (Should genes for which effect estimation fails be removed, marked with NAs or analysed by adding noise to them) 
# PARAMETER random.noise: "Random noise" TYPE DECIMAL FROM 0 TO 10000 DEFAULT 0.00001 (Random noise added to gene expression values. Use very small values like 0.00001 or less)

# JTT: 12.7.2006 Crated linear Mixed Model
# JTT: 19.10.2006: Modified to use nlme library
# MK: 22.05.2013: Noise factor parameter added 

#column.groups<-"group"
#column.random<-"age"

# Loads the libraries
library(nlme)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Test needs a parameter "groups" that specifies the grouping of the samples
# and a parameter random that specifies the random effect such as day or technician
groups<-phenodata[,which(column.groups==colnames(phenodata))]
random<-phenodata[,which(column.random==colnames(phenodata))]

# Sanity check
if(length(unique(random))==1) {
   stop("CHIPSTER-NOTE: You only have one level in variable random. Specify at least two levels for the random effect")
}

# Fits a linear mixed model for every gene, and saves residuals into a table
# Assumes no interaction between random and groups
# If estimation fails for a particular gene, return NAs
dat3<-dat2
error.vec <- rep(0, nrow(dat3))
for(i in 1:nrow(dat3)) {
	gene <- as.vector(as.numeric(dat3[i,]))
	residuals <- try(resid(lme(fixed=gene~groups, random=gene~1|random, control=list(maxIter=10000), method="REML")), silent=T)
	if(class(residuals) == "try-error") { 
		dat3[i,] <- rep(NA, length(gene));	
	} else {
		dat3[i,] <- residuals;
		error.vec[i] <- 1
	}
}

# Keep genes for which batch effects were estimated
if(error.handle == "remove") {
	dat3 <- dat3[which(error.vec==1),]
	calls <- calls[which(error.vec==1),]
}

# Re-analyse genes for which batch effects could not be estimated by adding random noise to their gene expression estimates 
if(error.handle == "add.noise") {
	for(i in which(error.vec == 0)) {
		residual.counter = 0
		class(residuals) <- "try-error"
		while(residual.counter < 1000 & class(residuals) == "try-error") {
			gene <- as.vector(as.numeric(dat2[i,]))
			gene <- gene + rnorm(length(gene), 0, random.noise)
			residuals <- try(resid(lme(fixed=gene~groups, random=gene~1|random, control=list(maxIter=10000), method="REML")), silent=T)
			residual.counter <- residual.counter +1
			print(paste(i, residual.counter))
		}
		
		if(class(residuals) == "try-error") { 
			print(i)
			stop("LME failed to estimate random effects. Please choose another error handling method or increse noise-level")
		} else {
			dat3[i,] <- residuals;	
		}
		#dat3[i,] <- resid(lme(fixed=gene~groups, random=gene~1|random, control=list(maxIter=10000), method="REML"))
	}
}

# Writes a table of results
write.table(data.frame(round(dat3, digits=2), calls), file="normalized-lme.tsv", sep="\t", row.names=T, col.names=T, quote=F)
