# TOOL stat-timeseries.R: "Time series" (Analyses of time series data. Finds periodically expressed genes. For the ICA method, a standard deviation threshold for significantly differentially expressed genes is needed. If there are more than one replicate per time point, this tool will not work. For maSigPro:phenodata:description must be unique,add column for time information,default column group is used as edesign.abiotic column replicate,add at least one column for each experimental condition as described in http://www.bioconductor.org/packages/2.7/bioc/vignettes/maSigPro/inst/doc/maSigPro-tutorial.pdf for edesign.abiotic;all parameters are default except vars="all")
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT timeseries.tsv: timeseries.tsv 
# OUTPUT profiles.pdf: profiles.pdf 
# PARAMETER column: column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the time to test)
# PARAMETER analysis.type: analysis.type TYPE [periodicity: periodicity, ica: ica, maSigPro: maSigPro] DEFAULT periodicity (Analysis type)
# PARAMETER p.value.threshold: p.value.threshold TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER p.value.adjustment.method: p.value.adjustment.method TYPE [yes: yes, no: no] DEFAULT yes (Apply Benjamimi-Hochberg correction?)
# PARAMETER SD.for.ICA: SD.for.ICA TYPE DECIMAL FROM 0 TO 10 DEFAULT 2.0 (Standard deviation for ICA)
# PARAMETER k.for.maSigPro: k.for.maSigPro TYPE DECIMAL FROM 0 TO 1000 DEFAULT 9 (maSigPro see.genes k=9)
# PARAMETER degree.for.maSigPro: degree.for.maSigPro TYPE DECIMAL FROM 0 TO 1000 DEFAULT 2 (maSigPro make.design.matrix degree=2)
# PARAMETER rsq.for.maSigPro: rsq.for.maSigPro TYPE DECIMAL FROM 0 TO 10 DEFAULT 0.7 (maSigPro get.siggenes rsq=0.7)
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Analysis methods for timeseries
# JTT 21.7.2006

# OH 10.10.2012
# added maSigPro as alternative time series analysis
# parameter added: k.for.maSigPro
# see documentation of maSigPro www.bioconductor.org/packages/2.7/bioc/vignettes/maSigPro/inst/doc/maSigPro-tutorial.pdf
#

# OH 17.01.2013
# added maSigPro as alternative time series analysis
# parameter added: degree.for.maSigPro
# see documentation of maSigPro www.bioconductor.org/packages/2.7/bioc/vignettes/maSigPro/inst/doc/maSigPro-tutorial.pdf
#

# OH 17.01.2013
# added maSigPro as alternative time series analysis
# parameter added: rsq.for.maSigPro
# see documentation of maSigPro www.bioconductor.org/packages/2.7/bioc/vignettes/maSigPro/inst/doc/maSigPro-tutorial.pdf
#


# Parameter settings (default) for testing purposes
#column<-"time"
#analysis.type<-"periodicity"
#p.value.threshold<-0.05
#p.value.adjustment.method<-"yes"
#SD.for.ICA<-2.0
#image.width<-600
#image.height<-600


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
times<-phenodata[,grep(column, colnames(phenodata))]

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# How many replicates there are per time point?
if(length(times)>length(levels(as.factor(times)))) {
   repl<-c()
   for(i in 1:length(levels(as.factor(times)))) {
      repl<-c(repl, sum(as.numeric(grep(times[1], times, value=T)==times[1])))
   }
} else {
   repl<-rep(1, length(times))
}

if(analysis=="periodicity") {
  # Load the libraries
  library(e1071)
  library(GeneCycle)
  #library(fastICA)

  # Making a longitudinal object
  dat3<-as.longitudinal(t(dat2), repeats=repl, time=times)

  # Replacing missing values
  dat4<-t(na.omit(t(impute(dat3))))

   f<-fisher.g.test(dat4)
   if(multiple.correction=="yes") {
      # Estimates the proportion of null p-values, uses BH-method
      p.adj<-fdrtool(f)$lfdr
      dev.off()
   } else {
      p.adj<-f
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
   dat3<-as.longitudinal(t(dat2), repeats=repl, time=times)

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
	library(e1071)
	library("maSigPro")

	adjust="none"
	if( multiple.correction=="yes" ) {
		adjust="BH"
	}

	samplecount=length(colnames(dat2))
	if( samplecount!=length(unique(phenodata["description"][[1]])) ) {
		stop("Provide a unique description for each sample in column description of phenodata.")
	}

	orig.colnames=colnames(dat2)
	colnames(dat2)=phenodata["description"][[1]]
	dat4=na.omit(impute(dat2))

	edesign=data.frame(Time=times)
	rownames(edesign)=phenodata["description"][[1]]
	edesign["Replicate"]=phenodata["group"][[1]]

	conditions=colnames(phenodata)[grep("sample|original_name|chiptype|group|description",colnames(phenodata),invert=TRUE)]
	conditions=conditions[grep(column,conditions,invert=TRUE)]
	if( length(conditions)==0 ) {
		stop("Provide at least one condition as a new column in phenodata.")
	}

	for( condition in conditions ) {
		edesign[condition]=phenodata[condition][[1]]
	}

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

	fit <- p.vector(dat4,design,Q=Q,MT.adjust=adjust,min.obs=min.obs)

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


