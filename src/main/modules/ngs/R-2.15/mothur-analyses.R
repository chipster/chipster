# TOOL mothur-analyses.R: "Statistical analysis for marker gene studies" (Compares the diversity or abundance between groups using several ANOVA-type of analyses. Statistical tests work only for datasets that contain 2-3 groups. Makes also an RDA ordination plot and rank abundance and rarefaction curves. Requires both a count table and a phenodata as inputs. Count table is simple tab-delimited text file where rows are samples and columns are taxa.)
# INPUT counttable.tsv: "Count table" TYPE GENERIC
# INPUT phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT rank-abundance_rarefaction_RDA.pdf: rank-abundance_rarefaction_RDA.pdf
# OUTPUT OPTIONAL stat-results.txt: stat-results.txt
# PARAMETER decostandm: "Method for standardizing species abundance values" TYPE [total: total, normalize: normalize, pa: pa, chi.square: chi.square, hellinger: hellinger, log: log] DEFAULT hellinger (Method for standardizing species abundance values before running the RDA and statistical analyses.)



# 28.8.2013

# Parameters for testing purposes
#setwd("C:\\Users\\tuimaja\\Desktop\\merged")
#decostandm<-"hellinger"

# Read the input data
dat<-read.table("counttable.tsv", header=T, sep="\t", row.names=1)
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
#phenodata$group<-sample(c(1,2), nrow(phenodata), replace=T)

# Some sanity checks
if(all(is.na(phenodata$group))) {
   stop("Specify the biological groups in the group column in the phenodata!")
}
if(any(!as.character(phenodata$group) %in% c("1", "2", "3"))) {
   stop("Biological groups in the group column of the phenodata need to be coded with 1, 2 and 3. Other characters are not allowed.")
}

# Rearrange data
dat2<-dat[phenodata$sample,]
tdat2<-t(dat2)

# Load the libraries
library(vegan)
library(rich)
library(BiodiversityR)
library(pegas)
library(labdsv)

# Open a report PDF
pdf("rank-abundance_rarefaction_RDA.pdf")

# Rarefaction curves
if(length(unique(phenodata$group))==1) {
   plot(specaccum(dat2, method="rarefaction"), ci.type="polygon", ci.col="lightblue", ci.lty=0, lwd=4, col="blue", main="Rarefaction curve(s)")
}
if(length(unique(phenodata$group))==2) {
   r1<-specaccum(dat2[phenodata$group==1,], method="rarefaction")
   r2<-specaccum(dat2[phenodata$group==2,], method="rarefaction")
   plot(r1, ci.type="polygon", ci.col="#ADD8E688", ci.lty=0, lwd=4, col="blue", ylim=c(0, max(r1$richness, r2$richness)+1), xlim=c(0,max(length(r1$richness), length(r2$richness))), main="Rarefaction curve(s)")
   plot(r2, ci.type="polygon", ci.col="#FFB6C188", ci.lty=0, lwd=4, col="red", add=TRUE)
}
if(length(unique(phenodata$group))==3) {
   r1<-specaccum(dat2[phenodata$group==1,], method="rarefaction")
   r2<-specaccum(dat2[phenodata$group==2,], method="rarefaction")
   r3<-specaccum(dat2[phenodata$group==3,], method="rarefaction")
   plot(r1, ci.type="polygon", ci.col="#ADD8E688", ci.lty=0, lwd=4, col="blue", ylim=c(0, max(r1$richness, r2$richness, r3$richness)+1), xlim=c(0,max(length(r1$richness), length(r2$richness), length(r3$richness))), main="Rarefaction curve(s)")
   plot(r2, ci.type="polygon", ci.col="#FFB6C188", ci.lty=0, lwd=4, col="red", add=TRUE)
   plot(r3, ci.type="polygon", ci.col="#90EE9088", ci.lty=0, lwd=4, col="green", add=TRUE)
}
legend(x="bottomright", lwd=4, lty=1, col=c("blue", "red", "green"), legend=c("Group 1","Group 2","Group 3"))


# Rank abundance curves
opar<-par()
par(mfrow=c(2,2))
rankabunplot(rankabundance(dat2), main="Rank - abundance plot, all samples")
ymax<-max(rankabundance(dat2))
if(length(unique(phenodata$group))==2) {
   rankabunplot(rankabundance(dat2[phenodata$group==1,]), main="Rank - abundance plot, group 1", ylim=c(0,ymax))
   rankabunplot(rankabundance(dat2[phenodata$group==2,]), main="Rank - abundance plot, group 2", ylim=c(0,ymax))
}
if(length(unique(phenodata$group))==3) {
   rankabunplot(rankabundance(dat2[phenodata$group==1,]), main="Rank - abundance plot, group 1", ylim=c(0,ymax))
   rankabunplot(rankabundance(dat2[phenodata$group==2,]), main="Rank - abundance plot, group 2", ylim=c(0,ymax))
   rankabunplot(rankabundance(dat2[phenodata$group==3,]), main="Rank - abundance plot, group 3", ylim=c(0,ymax))
}
par(opar)

# RDA plot
dat3<-decostand(dat2, method="hellinger")
group<-phenodata$group
if(length(unique(phenodata$group))==1) {
   rda1<-rda(dat3)
   l<-round(c(rda1$CA$eig/rda1$tot.chi*100)[1:2], digits=2)
   plot(rda1, type="n", scaling=3, main=paste("RDA -", decostandm),  xlab=paste("PCA1 - ", l[1], "%", sep=""), ylab=paste("PCA2 - ", l[2], "%", sep=""))
   #text(rda1, display="bp", lwd=4)
   points(rda1, display="species", pch="+", scaling=3, col="grey75")
   points(rda1, display="sites", scaling=3, pch=16, cex=2, col=ifelse(group==1, "#0066CC", "#CC0000"))
   text(rda1, display="sites", scaling=3, cex=0.5)
   sc<-scores(rda1, choises=c(1,2), display="sp", scaling=3)
}

if(length(unique(phenodata$group))>=2) {
   rda1<-rda(dat3~group, data.frame(group))
   l<-round(c(rda1$CCA$eig/rda1$tot.chi*100, rda1$CA$eig[1]/rda1$tot.chi*100), digits=2)
   plot(rda1, type="n", scaling=3, main=paste("RDA -", decostandm),  xlab=paste("RDA1 - ", l[1], "%", sep=""), ylab=paste("PC1 - ", l[2], "%", sep=""))
   text(rda1, display="bp", lwd=4)
   points(rda1, display="species", pch="+", scaling=3, col="grey75")
   points(rda1, display="sites", scaling=3, pch=16, cex=2, col=ifelse(group==1, "#0066CC", "#CC0000"))
   text(rda1, display="sites", scaling=3, cex=0.5)
   sc<-scores(rda1, choises=c(1,2), display="sp", scaling=3)
   legend(x="topright", bty="n", legend=paste("p-value for group = ", round(anova(rda1, by="margin")$'Pr(>F)'[1], digits=3), sep=""))
}

# Close the report PDF
dev.off()


# Statistical test will only be run, if there are two or three groups
if(length(unique(phenodata$group))>=2) {
   group<-phenodata$group

   res1<-contribdiv(t(dat2))

   dd<-dist(dat3)
   res2<-adonis(dd~group, permutations=9999)
   res3<-amova(dd~group, nperm=9999)
   res4<-anova(betadisper(dd, group))

   ind<-indval(dat2, group, numitr=9999)
   isa<-isamic(dat2, group)

   sink("stat-results.txt")
      cat("\n\n\n")
      cat("######################################\n")
      cat("### CONTRIBUTED DIVERSITY ANALYSIS ###\n")
      cat("######################################\n")
      cat("\n\n\n")
      print(res1)
      cat("\n\n\n")
      cat("###############################################################################\n")
      cat("### Permutational Multivariate Analysis of Variance Using Distance Matrices ###\n")
      cat("###############################################################################\n")
      cat("\n\n\n")
      print(res2)
      cat("\n\n\n")
      cat("######################################\n")
      cat("### Analysis of Molecular Variance ###\n")
      cat("######################################\n")
      cat("\n\n\n")
      print(res3)
      cat("\n\n\n")
      cat("##################################################################\n")
      cat("### Multivariate homogeneity of groups dispersions (variances) ###\n")
      cat("##################################################################\n")
      cat("\n\n\n")
      print(res4)
      cat("\n\n\n")
      cat("###################################################\n")
      cat("### Dufrene-Legendre Indicator Species Analysis ###\n")
      cat("###################################################\n")
      cat("\n\n\n")
      print(ind)
      cat("\n\n\n")
      cat("######################################################################\n")
      cat("### Indicator Species Analysis Minimizing Intermediate Occurrences ###\n")
      cat("######################################################################\n")
      cat("\n\n\n")
      print(isa)
   sink()
}
