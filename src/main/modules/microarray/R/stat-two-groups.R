# TOOL stat-two-groups.R: "Two groups tests" (Tests for comparing the mean gene expression of two groups. LPE only works, if the whole normalized data is used, i.e., the data should not be filtered. Other than empiricalBayes might be slow, if run on unfiltered data.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT two-sample.tsv: two-sample.tsv 
# PARAMETER column: "Column" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER OPTIONAL pairing: "Pairing" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing which samples form pairs. This option should be used if you have, for example, monitored your samples before and after treatment, have patient-matched data or you have expression data at multiple tissue sites from the same individuals, etc. LPE, F-test and fast-t-test do not support pairing information.)
# PARAMETER test: "Test" TYPE [empiricalBayes: "empirical Bayes", fast-t-test: "fast t-test", t-test: t-test, F-test: F-test, Mann-Whitney: Mann-Whitney, LPE: LPE, RankProd: RankProd] DEFAULT empiricalBayes (Test type)
# PARAMETER p.value.adjustment.method: "p-value adjustment method" TYPE [none: none, Bonferroni: Bonferroni, Holm: Holm, Hochberg: Hochberg, BH: BH, BY: BY] DEFAULT BH (Multiple testing correction method)
# PARAMETER p.value.threshold: "p-value threshold" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off \\(0.05\)for significant results)
# PARAMETER show.na: "Show NA" TYPE [yes: yes, no:no] DEFAULT no (include results where p-value is NA)

# Two-group parametric and non-parametric tests
# JTT 4.7.2006
# OH, 7.11.2011
# EK, 8.1.2012
# JT, 28.11.2012, fixed Wilcoxon test, sped up other tests
# MK, 05.02.2013, added paired tests for limma, t-test and Wilcox

# Loads the libraries
library(multtest)

# Renaming variables
meth<-test
if(test=="empiricalBayes" & (p.value.adjustment.method!="BH" & p.value.adjustment.method!="BY") ) {
	adj.method<-tolower(p.value.adjustment.method)
} else {
	adj.method<-p.value.adjustment.method
}
p.cut<-p.value.threshold
show.p.na<-show.na

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls <- dat[,grep("flag", names(dat))]
dat2  <- as.matrix(dat[,grep("chip", names(dat))])

# Test needs a parameter "groups" that specifies the grouping of the samples
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
groups<-phenodata[,pmatch(column,colnames(phenodata))]

if(exists("pairing")) {
	if(column == pairing) {
		stop("CHIPSTER-NOTE: Phenodata column describing the groups to test cannot be the same that contains pairing information")	
	}
	if(pairing!="EMPTY") {
		pairs<-phenodata[,pmatch(pairing,colnames(phenodata))]
	}
} else {
	pairing <- "EMPTY"
}

# Sanity checks
if(length(unique(groups))==1 | length(unique(groups))>=3) {
	stop("CHIPSTER-NOTE: You need to have exactly two groups to run this analysis")
}

# Testing

# Empirical Bayes
if(meth=="empiricalBayes") {
	library(limma)
	if(pairing=="EMPTY") {
		design<-model.matrix(~as.factor(groups))
	} else {
		groups		<- factor(groups)
		pairs		<- factor(pairs)
		#the 0 means no intercept for the first factor. In this case pairs-factor must be given first
		#design		<- model.matrix(~0+pairs+groups)
		design		<- model.matrix(~groups+pairs)
	}

	#stop(as.factor(groups));
	
	fit <- lmFit(dat2, design)
	fit <- eBayes(fit)
	tab <- toptable(fit, coef=2, number=nrow(fit), adjust.method=adj.method, sort.by="none")
	rows <- which(tab$adj.P.Val<=p.cut)
	p <- tab$adj.P.Val[rows]
	M <- tab$logFC[rows]
    dat <- data.frame(dat[rows,], p.adjusted=round(p, digits=6), FC=M);
	dat <- dat[order(dat$p.adjusted, decreasing=F),]

	write.table(dat, file="two-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

if(meth=="RankProd") {
	library(RankProd)
	
	dat2.1 <-dat2[,groups==unique(groups)[1]]
	dat2.2 <-dat2[,groups==unique(groups)[2]]
			
	if(pairing =="EMPTY") {
		group_vec <- c(rep(0, ncol(dat2.1)), rep(1, ncol(dat2.1)));
		dat.rp <- cbind(dat2.1, dat2.2);
		rp.fold.change <- apply(dat2.1, 1, mean, na.rm=T) - apply(dat2.2, 1, mean, na.rm=T)
	} else {
		pairs.1 <-pairs[groups==unique(groups)[1]]
		pairs.2 <-pairs[groups==unique(groups)[2]]
		
		#find shared elements (remove elements not present in one of the arrays)
		dat2.1 <- dat2.1[, which(pairs.1 %in% intersect(pairs.1, pairs.2))]
		dat2.2 <- dat2.2[, which(pairs.2 %in% intersect(pairs.1, pairs.2))]
		
		pairs.1 <- pairs.1[which(pairs.1 %in% intersect(pairs.1, pairs.2))];
		pairs.2 <- pairs.2[which(pairs.2 %in% intersect(pairs.1, pairs.2))];
		
		#average over those that are found multiple times in one
		temp.1 <- t(rowsum(t(dat2.1), factor(pairs.1), reorder = FALSE, na.rm = TRUE));
		temp.2 <- t(rowsum(1L - is.na(t(dat2.1)), factor(pairs.1), reorder = FALSE));
		dat2.1 <- temp.1 / temp.2;
		
		temp.1 <- t(rowsum(t(dat2.2), factor(pairs.2), reorder = FALSE, na.rm = TRUE));
		temp.2 <- t(rowsum(1L - is.na(t(dat2.2)), factor(pairs.2), reorder = FALSE));
		dat2.2 <- temp.1 / temp.2;

		#sort columns so that samples have the same order
		dat2.2 <- dat2.2[,match(pairs.1, pairs.2)]
		
		if(ncol(dat2.1) != ncol(dat2.2)) { stop("Paired RankProd error: matrices differn in column number")}
		dat.rp <- dat2.1 - dat2.2;
		group_vec <- rep(0, ncol(dat.rp));

		rp.fold.change <- apply(dat.rp, 1, mean, na.rm=T)
	}

	RPdata 	<- RP(dat.rp, cl=group_vec, num.perm=10, logged=TRUE)
	p.raw   <- cbind(RPdata$pval[,1],  RPdata$pval[,2]);
}
		
# Fast T-test
if(meth=="fast-t-test") {
	if(pairing!="EMPTY") { stop("CHIPSTER-NOTE: Fast T-test does not support pairing information"); }
		
	fit1<-lm(t(dat2)~groups)
	p<-rep(NA, nrow(dat2))
	for(i in 1:nrow(dat2)) {
		sum(fit1$residuals[,i]^2)->sse
		sum((dat2[i,]-mean(as.numeric(dat2[i,])))^2)->sst
		r2<-1-(sse/sst)
		f<-r2/((1-r2)/(ncol(dat2)-1))
		p[i]<-1-pf(f, 1, (ncol(dat2)-1))
	}
	p.raw<-p	
}

# T-test
if(meth=="t-test") {
	dat2.1 <-dat2[,groups==unique(groups)[1]]
	dat2.2 <-dat2[,groups==unique(groups)[2]]
	p      <- rep(as.numeric(NA), nrow(dat2))
	if (pairing=="EMPTY") {
		for(i in 1:nrow(dat2)) {
			if( (sum(!is.na(dat2.1[i,])) > 1) & (sum(!is.na(dat2.2[i,])) > 1) &&
				( length(which(!dat2.1[i,]==mean(dat2.1[i,])))>0 || length(which(!dat2.2[i,]==mean(dat2.2[i,])))>0 )
			  ) { 
				p[i] <- t.test(x=dat2.1[i,], y=dat2.2[i,])$p.value
			} 
		}
	} else {
		pairs.1 <-pairs[groups==unique(groups)[1]]
		pairs.2 <-pairs[groups==unique(groups)[2]]

		#find shared elements (remove elements not present in one of the arrays)
		dat2.1 <- dat2.1[, which(pairs.1 %in% intersect(pairs.1, pairs.2))]
		dat2.2 <- dat2.2[, which(pairs.2 %in% intersect(pairs.1, pairs.2))]

		pairs.1 <- pairs.1[which(pairs.1 %in% intersect(pairs.1, pairs.2))];
		pairs.2 <- pairs.2[which(pairs.2 %in% intersect(pairs.1, pairs.2))];
	
		#average over those that are found multiple times in one
		temp.1 <- t(rowsum(t(dat2.1), factor(pairs.1), reorder = FALSE, na.rm = TRUE));
		temp.2 <- t(rowsum(1L - is.na(t(dat2.1)), factor(pairs.1), reorder = FALSE));
		dat2.1 <- temp.1 / temp.2;

		temp.1 <- t(rowsum(t(dat2.2), factor(pairs.2), reorder = FALSE, na.rm = TRUE));
		temp.2 <- t(rowsum(1L - is.na(t(dat2.2)), factor(pairs.2), reorder = FALSE));
		dat2.2 <- temp.1 / temp.2;

		#sort columns so that samples have the same order
		dat2.2 <- dat2.2[,match(pairs.1, pairs.2)]

		for(i in 1:nrow(dat2)) {
			print(dat2.1[i,])
			print(dat2.2[i,])
			
			if((sum(!is.na(dat2.1[i,])) > 1) & (sum(!is.na(dat2.2[i,])) > 1) && length(unique(round(dat2.1[i,] - dat2.2[i,], digits=10))) > 1 &&
				( length(which(!dat2.1[i,]==mean(dat2.1[i,])))>0 || length(which(!dat2.2[i,]==mean(dat2.2[i,])))>0 )
		    ) { 
				print(dat2.1[i,])
				print(dat2.2[i,])
				p[i] <- t.test(x=dat2.1[i,], y=dat2.2[i,], paired=TRUE)$p.value
			}
		}	
	}

	p.raw<-p
}

# F-test
if(meth=="F-test") {
	#paired var.test does in fact exist in package PairedData
	if(pairing!="EMPTY") { stop("CHIPSTER-NOTE: F-test does not support pairing information"); }
	
	dat2.1 <-dat2[,groups==unique(groups)[1]]
	dat2.2 <-dat2[,groups==unique(groups)[2]]
	p 	   <- rep(as.numeric(NA), nrow(dat2))
	for(i in 1:nrow(dat2)) {
		p[i] <- var.test(x=as.numeric(dat2.1[i,]), y=as.numeric(dat2.2[i,]))$p.value
	}
	p.raw<-p
	
}

# Mann-Whitney test
if(meth=="Mann-Whitney") {
	dat2.1 <-dat2[,groups==unique(groups)[1]]
	dat2.2 <-dat2[,groups==unique(groups)[2]]
	p      <- rep(as.numeric(NA), nrow(dat2))

	if (pairing=="EMPTY") {
		for(i in 1:nrow(dat2)) {
			p[i] <- wilcox.test(x=dat2.1[i,], y=dat2.2[i,])$p.value
		}
	} else {
		pairs.1 <-pairs[groups==unique(groups)[1]]
		pairs.2 <-pairs[groups==unique(groups)[2]]
		
		#find shared elements (remove elements not present in one of the arrays)
		dat2.1 <- dat2.1[, which(pairs.1 %in% intersect(pairs.1, pairs.2))]
		dat2.2 <- dat2.2[, which(pairs.2 %in% intersect(pairs.1, pairs.2))]
		
		pairs.1 <- pairs.1[which(pairs.1 %in% intersect(pairs.1, pairs.2))];
		pairs.2 <- pairs.2[which(pairs.2 %in% intersect(pairs.1, pairs.2))];
		
		#average over those that are found multiple times in one
		temp.1 <- t(rowsum(t(dat2.1), factor(pairs.1), reorder = FALSE, na.rm = TRUE));
		temp.2 <- t(rowsum(1L - is.na(t(dat2.1)), factor(pairs.1), reorder = FALSE));
		dat2.1 <- temp.1 / temp.2;
		
		temp.1 <- t(rowsum(t(dat2.2), factor(pairs.2), reorder = FALSE, na.rm = TRUE));
		temp.2 <- t(rowsum(1L - is.na(t(dat2.2)), factor(pairs.2), reorder = FALSE));
		dat2.2 <- temp.1 / temp.2;

		#sort columns so that samples have the same order
		dat2.2 <- dat2.2[,match(pairs.1, pairs.2)]
		
		for(i in 1:nrow(dat2)) {
			p[i] <- wilcox.test(x=dat2.1[i,], y=dat2.2[i,], paired=TRUE)$p.value
		}
	}
		
	p.raw<-p
}

## these methods have a common p-value adjustment and result generation
if (meth %in% c("Mann-Whitney", "fast-t-test", "t-test", "F-test")){
	if(!exists("p.raw")) {	
		stop("Could not find p.raw object")
	}	
	
	if (adj.method %in% c("Bonferroni", "Holm", "Hochberg")) {
		p.adjusted <- p.adjust(p.raw, method=tolower(adj.method))
	} else {
		p.adjusted <- p.adjust(p.raw, method=adj.method)
	}
	if( show.p.na == "yes" ) {
		dat<-dat[which(p.adjusted<=p.cut|is.na(p.adjusted)),]   
		p.adjusted<-p.adjusted[which(p.adjusted<=p.cut|is.na(p.adjusted))]
	} else {
		dat<-dat[which(p.adjusted<=p.cut),]   
		p.adjusted<-p.adjusted[which(p.adjusted<=p.cut)]
	}
	write.table(data.frame(dat, p.adjusted=round(p.adjusted, digits=6)),
			file="two-sample.tsv", sep="\t", row.names=TRUE,
			col.names=TRUE, quote=FALSE)
}

if(meth=="RankProd") {
	if(!exists("p.raw")) {	
		stop("Could not find p.raw object")
	}	
	
	p.adjusted <- matrix(NA, ncol=ncol(p.raw), nrow=nrow(p.raw));
	if (adj.method %in% c("Bonferroni", "Holm", "Hochberg")) {
		p.adjusted[, 1] <- p.adjust(p.raw[, 1], method=tolower(adj.method));
		p.adjusted[, 2] <- p.adjust(p.raw[, 2], method=tolower(adj.method));
	} else {
		p.adjusted[, 1] <- p.adjust(p.raw[, 1], method=adj.method)
		p.adjusted[, 2] <- p.adjust(p.raw[, 2], method=adj.method)
	}
	if( show.p.na == "yes" ) {
		dat <- dat[which(apply(p.adjusted,1,min)<=p.cut|is.na(apply(p.adjusted,1,min))),]   
		p.adjusted <- p.adjusted[which(apply(p.adjusted,1,min)<=p.cut|is.na(apply(p.adjusted,1,min))), ]
		p.adjusted <- apply(p.adjusted,1,min)

		rp.fold.change<-rp.fold.change[which(p.adjusted<=p.cut|is.na(p.adjusted))]  
	} else {
		dat <- dat[which(apply(p.adjusted,1,min)<=p.cut),]   
		p.adjusted <- p.adjusted[which(apply(p.adjusted,1,min)<=p.cut), ]
		p.adjusted <- apply(p.adjusted,1,min)

		rp.fold.change<-rp.fold.change[which(p.adjusted<=p.cut)]
	}
	
	write.table(data.frame(dat, p.adjusted=round(p.adjusted, digits=6), FC=round(rp.fold.change, digits=2)),
			file="two-sample.tsv", sep="\t", row.names=TRUE,
			col.names=TRUE, quote=FALSE)
}
	
if(meth=="LPE") {
	if(pairing!="EMPTY") { stop("CHIPSTER-NOTE: LPE does not support pairing information"); }
	
	library(LPE)

	group1 <-dat2[,groups==unique(groups)[1]]
	group2 <-dat2[,groups==unique(groups)[2]]	
	g1.x<-baseOlig.error(group1)
	g2.x<-baseOlig.error(group2)
	lp<-data.frame(lpe(group1, group2, g1.x, g2.x, probe.set.name=row.names(dat2)))

	if(adj.method=="none" | adj.method=="Holm") {
		x.location <- grep("^x", names(lp))
		y.location <- grep("^y", names(lp))
		x <- lp[, x.location]
		y <- lp[, y.location]
		pnorm.diff <- pnorm(lp$median.diff, mean = 0, sd = lp$pooled.std.dev)
		
		p.adjusted <- 2 * apply(cbind(pnorm.diff, 1 - pnorm.diff), 1, min)
		
		if(adj.method=="Holm") {
			#do Holm-Bonferroni correction
			p.adjusted <- pmin(p.adjusted * seq(1,length(p.adjusted)),1) 	
		}
	}
	if(adj.method=="Hochberg") {
		stop("CHIPSTER-NOTE: LPE function does not support Hochberg's correction. Please rerun the function again after choosing another p.value.adjustment.method")
	}
	if(adj.method=="Bonferroni" | adj.method=="BH" | adj.method=="BY") {
		p.adjusted <- fdr.adjust(lp, adjp=adj.method)
		p.adjusted <- p.adjusted[order(match(rownames(p.adjusted),rownames(lp))),]
		p.adjusted <- data.frame(p.adjusted)$FDR;
		
		#dat2<-merge(dat,as.data.frame(round(fdr, digits=4)), by.x="row.names", by.y="row.names")
		#dat2<-dat2[,-ncol(dat2)] # Removes the last columns that holds Z-test values
		#names(dat2)[which(names(dat2)=="FDR")]<-"p.adjusted" # Renames "FDR" with "p.adjusted"
	}
	if( show.p.na == "yes" ) {
		dat<-dat[which(p.adjusted<=p.cut|is.na(p.adjusted)),]   
		p.adjusted<-p.adjusted[which(p.adjusted<=p.cut|is.na(p.adjusted))]
	} else {
		dat<-dat[which(p.adjusted<=p.cut),]   
		p.adjusted<-p.adjusted[which(p.adjusted<=p.cut)]
	}
	write.table(data.frame(dat, p.adjusted=round(p.adjusted, digits=6)), file="two-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

# EOF
