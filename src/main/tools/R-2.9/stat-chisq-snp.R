# ANALYSIS Statistics/"Association analysis" (Association tests for normalized SNP array data. Runs a Chi square test
# for every SNP. Hardy-Weinberg equilibrium is tested in controls only. Association tests use the grouping information
# of sample in group column of phenodata. Association tests are run for genotype frequences and dominant and recessive
# models.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT assoc.tsv
# PARAMETER test.for [Hardy-Weinberg, association] DEFAULT association (What to test for)
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to compare)


# Association analysis with normalized SNP data
# 25.4.2008
#
# modified by MG, 14.10.2009

# Read in data
dat<-read.table("normalized.tsv", header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Test needs a parameter "groups" that specifies the grouping of the samples
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
groups<-phenodata[,pmatch(column,colnames(phenodata))]

# Select only samples in group 1 (controls)
dat2<-dat[,which(groups==1)]

# Testing for Hardy-Weinberg in controls
if(test.for=="Hardy-Weinberg") {
   p<-rep("NA", nrow(dat))
   for(i in 1:nrow(dat)) {
      table1<-table(as.numeric(as.vector(dat2[i,])))
      if(length(table1)==1) {
         p[i]<-1
      }
      if(length(table1)==2) {
        p[i]<-1
      }
      if(length(table1)==3) {
         AA<-table1[1]
         AB<-table1[2]
         BB<-table1[3]
         s<-sum(table1)
         fa<-(AA*2+AB)/(2*s)
         fb<-(BB*2+AB)/(2*s)  
         o<-c(AA,AB,BB)
         e<-c(fa^2*s, 2*fa*fb*s, fb^2*s)
         p[i]<-chisq.test(x=o, y=e, simulate.p.value=T)$p.value
      }
   }
  # Saving the p-values
  p.hardy<-p
}

# Excluding the SNPs that violate Hardy-Weinberg
# dat2<-dat[,grep("chip", names(dat))]
# dat2<-dat2[p>=0.05,]
# dat2<-as.matrix(dat2)

# Testing for association
if(test.for=="association") {
   # Testing for difference between two groups
   library(scrime)
   
   
############################
# MG ADDED THESE TWO LINES #
############################

   dat2<-dat[,grep("chip", names(dat))]
   dat2<-as.matrix(dat2)
   
############################
   
   
   
   p<-rowChisqStats(dat2, cl=groups, compPval = TRUE)
   p<-p$rawp

   # Saving the p-values
   p.genotype<-p

   # Getting the raw data again
   dat2<-dat[,grep("chip", names(dat))]
   dat2<-as.matrix(dat2)

   # Testing the dominant model
   dat2[dat2==2]<-1
   dat2[dat2==3]<-2
   p<-rowChisqStats(dat2, cl=groups, compPval = TRUE)
   p<-p$rawp

   # Saving the p-values
   p.dominant<-p
  
   # Getting the raw data again
   dat2<-dat[,grep("chip", names(dat))]
   dat2<-as.matrix(dat2)

   # Testing the recessive model
   dat2[dat2==2]<-3
   dat2[dat2==3]<-2
   p<-rowChisqStats(dat2, cl=groups, compPval = TRUE)
   p<-p$rawp
   
   # Saving the p-values
   p.recessive<-p
}

# Writing data to disk
if(test.for=="Hardy-Weinberg") {
   write.table(data.frame(dat, p.adjusted=p.hardy), file="assoc.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}
if(test.for=="association") {
   write.table(data.frame(dat, p.genotype, p.dominant, p.recessive), file="assoc.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}


