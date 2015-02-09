# TOOL stat-timeseries.R: "Time series" (Analyses of time series data. The ICA and periodicity methods can handle only single time course microarray experiments, whereas maSigPro is suitable also for the analysis of multiseries time course microarray experiments.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT timeseries.tsv: timeseries.tsv 
# OUTPUT profiles.pdf: profiles.pdf 
# PARAMETER analysis.type: "Analysis type" TYPE [periodicity: periodicity, ica: ICA, maSigPro: maSigPro] DEFAULT periodicity (Analysis type)
# PARAMETER column: "Time column" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the time to test)
# PARAMETER rep.column: "Replicate column for maSigPro" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column indicating the replicate groups)
# PARAMETER other.start: "First experimental column for maSigPro" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata should include columns that give the assignment of arrays to experimental groups. There should be as many columns as experimental groups and these columns should reside side by side. Set this parameter to the first such column)
# PARAMETER other.end: "Last experimental column for maSigPro" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata should include columns that give the assignment of arrays to experimental groups. Set this column to the last experimental columns)
# PARAMETER k.for.maSigPro: "k for maSigPro" TYPE DECIMAL FROM 0 TO 1000 DEFAULT 9 (maSigPro see.genes k=9)
# PARAMETER degree.for.maSigPro: "Degree for maSigPro" TYPE DECIMAL FROM 0 TO 1000 DEFAULT 2 (maSigPro make.design.matrix degree=2)
# PARAMETER rsq.for.maSigPro: "rsq for maSigPro" TYPE DECIMAL FROM 0 TO 10 DEFAULT 0.7 (maSigPro get.siggenes rsq=0.7)
# PARAMETER SD.for.ICA: "SD for ICA" TYPE DECIMAL FROM 0 TO 10 DEFAULT 2.0 (Standard deviation for ICA)
# PARAMETER OPTIONAL p.value.threshold: "p-value threshold" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER OPTIONAL p.value.adjustment.method: "p-value adjustment method" TYPE [none: none, Bonferroni: Bonferroni, Holm: Holm, Hochberg: Hochberg, BH: BH, BY: BY] DEFAULT BH (Multiple testing correction method)
# PARAMETER OPTIONAL image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER OPTIONAL image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# JTT 21.7.2006: Analysis methods for timeseries
# OH 10.10.2012: added maSigPro as alternative time series analysis, parameter added: k.for.maSigPro
# OH 17.01.2013: added maSigPro as alternative time series analysis, parameter added: degree.for.maSigPro
# OH 17.01.2013: added maSigPro as alternative time series analysis, parameter added: rsq.for.maSigPro
# MK 15.11.2013: bugs fixed in ICA and periodicity and maSigPro interface polished up. Manual upadated for maSigPro

# Renaming variables
analysis<-analysis.type
p.cut<-p.value.threshold
multiple.correction<-p.value.adjustment.method
thresh<-SD.for.ICA
w<-image.width 
h<-image.height 

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)
symbols=list()
if( length(grep("symbol",colnames(dat)))>0 ) {
	symbols=dat["symbol"]
}
description=list()
if( length(grep("description",colnames(dat)))>0 ) {
	description=dat["description"]
}

# Reads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
times<-phenodata[column][[1]]

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# How many replicates there are per time point?
repl <- as.vector(table(times))[match(unique(times), dimnames(table(times))$times)]

if(analysis=="periodicity") {
	# Load the libraries
	library(e1071)
	library(GeneCycle)
	
	# Making a longitudinal object
	dat3<-as.longitudinal(t(dat2), repeats=repl, time=unique(times))
	
	# Replacing missing values
	dat4<-t(na.omit(t(impute(dat3))))
	
	f<-fisher.g.test(dat4)
	
	if (p.value.adjustment.method %in% c("Bonferroni", "Holm", "Hochberg")) {
		p.adj <- p.adjust(f, method=tolower(p.value.adjustment.method))
	} else {
		p.adj <- p.adjust(f, method=p.value.adjustment.method)
	}
	
	dat5<-as.data.frame(t(dat4))
	names(dat5)<-names(dat2)
	write.table(data.frame(dat5[p.adj<=p.cut,], p.adjusted=round(p.adj[p.adj<=p.cut], digits=6)), file="timeseries.tsv", sep="\t", row.names=T, col.names=T, quote=F)
	pdf(file="profiles.pdf", width=w/72, height=h/72)
	plot(1, 1, col=0)
	text(1, 1, "This is a dummy image.", col=1)
	text(1, 0.9, "To generate an image of gene expression profiles, use ica option.", col=1)
	dev.off()
}

if(analysis=="ica") {
	# Load the libraries
	library(e1071)
	library(GeneCycle)
	library(fastICA)
	
	# Making a longitudinal object
	dat3<-as.longitudinal(t(dat2), repeats=repl, time=unique(times))
	
	# Replacing missing values
	dat4<-t(na.omit(t(impute(dat3))))
	
	# Calculating independent component analysis usign a fast method
	o<-fastICA(t(dat4), n.comp=(ncol(t(dat4))-1), method="C")
	d<-c()
	g<-c()
	for(i in 1:ncol(o$S)) {
		s<-o$S[,i]
		l<-which(s>=(mean(s)+thresh*sd(s)))
		d<-c(d, l)
		g<-c(g, rep(i, length(l)))
		dg<-data.frame(dat[d,], cluster=g)
	}
	write.table(dg, file="timeseries.tsv", sep="\t", row.names=T, col.names=T, quote=F)
	a<-data.frame(times, t(o$A))
	pdf(file="profiles.pdf", width=w/72, height=h/72)
	par(mar=c(0, 1, 1, 0)+0.1)
	par(mfrow=c(ceiling(sqrt(ncol(a))), ceiling(sqrt(ncol(a)))))
	for(i in 2:ncol(a)) {
		plot(a$time, a[,i], xaxt="n", type="l", xlab=NULL, ylab=NULL, main=paste("Chip", i-1, sep=" "), cex.main=0.75, cex.axis=0.7, cex.lab=0.75, tck=-0.01, mgp=c(3,0.2,0))
	}
	dev.off()
}

if(analysis=="maSigPro") {
	suppressPackageStartupMessages(suppressWarnings(library(tcltk)))
	library(e1071)
	library("maSigPro")
	
	adjust <- p.value.adjustment.method
	
	orig.colnames=colnames(dat2)
	colnames(dat2) <- make.names(phenodata["description"][[1]], unique=T)
	dat4=na.omit(impute(dat2))
	
	edesign=data.frame(Time=times)
	edesign["Replicate"]=phenodata[rep.column][[1]]
	edesign <- cbind(edesign, phenodata[grep(other.start, colnames(phenodata)):grep(other.end, colnames(phenodata))])
	rownames(edesign)=colnames(dat2)
	
	Q=p.cut
	alfa=Q
	cluster.data = 1
	k = k.for.maSigPro
	degree = degree.for.maSigPro
	rsq=rsq.for.maSigPro
	cluster.method = "hclust"
	distance = "cor"
	agglo.method = "ward"
	summary.mode = "median"
	color.mode = "rainbow"
	show.fit = TRUE
	step.method = "backward"
	nvar.correction = FALSE
	min.obs = 3
	show.lines = TRUE
	cexlab = 0.8
	legend = TRUE
	#vars="groups"
	vars="all"
	main=NULL
	
	design <- make.design.matrix(edesign, degree = degree)
	dis <- design$dis
	
	fit <- p.vector(dat4, design, Q=Q, MT.adjust=adjust, min.obs=min.obs)
	if(min(fit$p.adjust) > Q) {
		assign("last.warning", NULL, envir = baseenv())
		stop("CHIPSTER-NOTE: No significant genes found.")	
	}
	
	tstep <- T.fit(fit,step.method=step.method,alfa=alfa)
	
	sigs <- get.siggenes(tstep,rsq=rsq,vars=vars)
	
	summary <- sigs$summary
	sig.genes <- sigs$sig.genes
	sig.genes <- sig.genes
	
	#from function maSigPro
	pdf(file = "profiles.pdf")
	
	if (!is.null(sig.genes)) {
		#for (i in 1:length(sig.genes)) {
		if (nrow(sig.genes[[1]]) > 0) {
			print("running see.genes")
			cluster <- see.genes(data = sig.genes,
					cluster.data = cluster.data, k = k, cluster.method = cluster.method,
					distance = distance, agglo.method = agglo.method,
					show.fit = show.fit, dis = dis, step.method = step.method,
					min.obs = min.obs, alfa = alfa, nvar.correction = nvar.correction,
					summary.mode = summary.mode, color.mode = color.mode,
					show.lines = show.lines, cexlab = cexlab,
					legend = legend, newX11 = FALSE, main = main)
			sig.genes[[1]] <- cbind(sig.genes[[1]],cluster$cut)
			cluster.algorithm <- cluster$cluster.algorithm.used
			groups <- cluster$groups
		}
		#}
	}
	
	dev.off()
	
	output <- list(summary, sig.genes, fit$dat, fit$G, edesign, dis, fit$min.obs, fit$p.vector, tstep$variables, tstep$g, fit$BH.alfa, step.method, Q, alfa, tstep$influ.info, vars, cluster.algorithm, groups)
	names(output) <- c("summary", "sig.genes", "input.data", "G", "edesign", "dis", "min.obs", "p.vector", "variables", "g", "BH.alfa", "step.method", "Q", "step.alfa", "influ.info", "select.vars", "cluster.algorithm.used", "groups")
	
	result = output$sig.genes$sig.profiles
	colnames(result)[1:length(orig.colnames)]=orig.colnames
	result["p-value"]=round(output$sig.genes$sig.pvalues$"p-value",digits=6)
	
	if( length(description)>0 ) {
		result2=merge(description,result,by = "row.names")
		rownames(result2)=result2[[1]]
		result=result2[-1]
	}
	if( length(symbols)>0 ) {
		result2=merge(symbols,result,by = "row.names")
		rownames(result2)=result2[[1]]
		result=result2[-1]
	}
	
	write.table(result,file="timeseries.tsv",sep="\t",row.names=T,col.names=T,quote=F)
	
}


