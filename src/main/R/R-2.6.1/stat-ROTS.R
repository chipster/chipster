# ANALYSIS Statistics/ROTS (Reproducibility-optimized test statistic for comparing the expression levels between two groups; Elo et al. 2008) 
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT ROTS.tsv
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER fdr.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (FDR cut-off for significant results)
# PARAMETER B INTEGER FROM 0 TO 100000 DEFAULT 500 (Number of resamplings)
# PARAMETER K INTEGER FROM 1000 TO 100000 DEFAULT 5000 (Largest top list size considered)



# Reproducibility-optimized two-group test

# LLE 11.9.2008



# Loads the library needed for fitting a linear model to the expression data matrix.

library(SAGx)



# Loads the normalized data

file <- c("normalized.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)



# Separates expression values and flags

calls <- dat[,grep("flag", names(dat))]
data <- dat[,grep("chip", names(dat))]



# Test needs a parameter "groups" that specifies the grouping of the samples

phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
groups <- phenodata[,grep(column, colnames(phenodata))]



# Sanity checks

if(length(unique(groups))==1 | length(unique(groups))>=3) {
   stop("You need to have exactly two groups to run this analysis")
}



# Bootstrap samples

bootstrapSamples<- function(B, labels) {

   l1 <- which(labels==1)
   l2 <- which(labels==2)

   samples <- matrix(nrow=B, ncol=length(l1)+length(l2))
   for(i in 1:B) {

      j <- 1
      while(j==1) {
         a1 <- sample(l1, length(l1), replace=TRUE)
	 j <- length(unique(a1))
      }

      j <- 1
      while(j==1) {
         a2 <- sample(l2, length(l2), replace=TRUE)
	 j <- length(unique(a2))
      }

      samples[i,] <- c(a1,a2)
   }

   return(samples)	
}	



# Overlap calculations (N = vector of top list sizes)

calculateOverlap<- function(result1, result2, N) {
   
   temp <- abs(result2)[ order( abs(result1), decreasing = T )]
   temp2 <- sort( abs(result2), decreasing = T )
   top<- sapply( N, function(i){ sum( temp[1:i] >= temp2[i] ) } )
   
   return(top / N)
}



# FDR calculations

biggerEq <- function(x, y) {
 
   x <- sort(x, decreasing = TRUE)		
   y <- sort(y, decreasing = TRUE)		
   a <- match(x, x)				
   b <- x %in% y				
   z <- sort(c(x, y), decreasing = TRUE)		
   d <- match(x, z)				
 
   return(d - a + b)
}


calculateFDR<- function(observed, permuted) {
   
   observed <- abs(observed)
   permuted <- abs(permuted)
   ord <- order(observed, decreasing = TRUE)
   a <- observed[ord]
   
   A <- matrix(NA, nrow=length(a), ncol=ncol(permuted))
   for(i in 1:ncol(A)) {
      a.rand <- sort(permuted[,i], decreasing = TRUE)
      n.bigger <- biggerEq(a, a.rand)
      A[ord,i] <- n.bigger/(1:length(a))
   }

   FDR <- apply(A, 1, median)
   FDR[FDR>1] <- 1
	
   return(FDR)
}



# Reproducibility-optimized test statistic

ROTS <- function(data, groups, B, K) {
   
   ssq <- c(c(0:20) / 100, c(11:50)/50, c(6:25) / 5)
   N <- c( c(1:20)*5, c(11:50)*10, c(21:40)*25, c(11:1000)*100)
   
   # The top list size cannot be larger than the total number of genes
   K <- min(K,nrow(data))					
   N <- N[N < K] 					
    
   # Reorder the data according to the group labels
   data1 <- data[,groups==unique(groups)[1]]
   data2 <- data[,groups==unique(groups)[2]]
   data <- cbind(data1,data2)
   cl <- c(rep(1, ncol(data1)), rep(2, ncol(data2)))
   rm(data1)
   rm(data2)
    
   samples <- bootstrapSamples(2*B, cl)

   design <- model.matrix(~as.factor(cl))
   contrast <- c(0, 1)

   D<- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
   S<- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
   a<- design
   for (i in 1:nrow(samples)) {
      samples.R <- samples[i,]
      a[, contrast > 0] <- design[samples.R, contrast > 0]
      fit <- Xprep(indata=data[, samples.R], formula=NULL, design=a, contrast=contrast)
      D[, i] <- fit$Mbar
      S[, i] <- sqrt(fit$Vest / fit$k)
   }

   pD<- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
   pS<- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
   a<- design
   for (i in 1:nrow(samples)) {
      samples.R <- sample(1:ncol(data))
      fit <- Xprep(indata=data[, samples.R], formula=NULL, design=a, contrast=contrast)
      pD[, i] <- fit$Mbar
      pS[, i] <- sqrt(fit$Vest / fit$k)
   }

   reprotable <- matrix(nrow=length(ssq) + 1, ncol=length(N))
   colnames(reprotable) <- N
   row.names(reprotable) <- c(ssq, "slr")

   reprotable.P <- matrix(nrow=length(ssq) + 1, ncol=length(N))
   colnames(reprotable.P) <- N
   row.names(reprotable.P) <- c(ssq, "slr")

   reprotable.sd <- matrix(nrow=length(ssq) + 1, ncol=length(N))
   colnames(reprotable.sd) <- N
   row.names(reprotable.sd) <- c(ssq, "slr")

   for(i in 1:length(ssq)) {
   
      overlaps <- matrix(nrow=B, ncol=length(N))
      overlaps.P <- matrix(nrow=B, ncol=length(N))
      for(b in 1:B) {
         result1<- D[, b]/(S[, b] + ssq[i])
	 result2<- D[, b + B]/(S[, b + B] + ssq[i])
         overlaps[b, ] <- calculateOverlap(result1, result2, N)
         result1<- pD[, b]/(pS[, b] + ssq[i])
         result2<- pD[, b + B]/(pS[, b + B] + ssq[i])
         overlaps.P[b, ] <- calculateOverlap(result1, result2, N)
      }
      reprotable[i, ] <- apply(overlaps, 2, mean)
      reprotable.P[i, ] <- apply(overlaps.P, 2, mean)
      reprotable.sd[i, ] <- apply(overlaps, 2, sd)

   }
   
   i <- length(ssq) + 1

      overlaps <- matrix(nrow=B, ncol=length(N))
      overlaps.P <- matrix(nrow=B, ncol=length(N))
      for(b in 1:B) {
         result1<- D[, b]
	 result2<- D[, b + B]
         overlaps[b, ] <- calculateOverlap(result1, result2, N)
         result1<- pD[, b]
         result2<- pD[, b + B]
         overlaps.P[b, ] <- calculateOverlap(result1, result2, N)
      }
      reprotable[i, ] <- apply(overlaps, 2, mean)
      reprotable.P[i, ] <- apply(overlaps.P, 2, mean)
      reprotable.sd[i, ] <- apply(overlaps, 2, sd)

   ztable <- (reprotable - reprotable.P) / reprotable.sd
   sel <- which(ztable == max(ztable), arr.ind=TRUE)
   if(length(sel)>2) sel<- sel[1,]
   if(sel[1] < nrow(reprotable)) {
      a1 <- as.numeric(row.names(reprotable)[sel[1]])
      a2 <- 1
   }
   if(sel[1] == nrow(reprotable)) {
      a1 <- 1
      a2 <- 0
   }
   k <- as.numeric(colnames(reprotable)[sel[2]])
   R <- reprotable[sel[1],sel[2]]
   Z <- ztable[sel[1],sel[2]]

   fit <- Xprep(indata=data, formula=~as.factor(cl), contrast=contrast)
   d <- t(fit$Mbar / (a1 + a2 * sqrt(fit$Vest / fit$k)))		
   FDR<- calculateFDR(d, pD/(a1 + a2 * pS))

   return(list(
      d=d,		
      FDR=FDR,		
      a1=a1,		
      a2=a2,		
      k=k,		
      R=R,		
      Z=Z		
   ))
    
}



# Testing

  result <- ROTS(data, groups, B, K) 
  d <- result$d
  FDR <- result$FDR
  FDR[FDR > fdr.threshold] <- NA
  write.table(na.omit(data.frame(dat, d=round(d, digits=2), FDR=round(FDR, digits=4))), file="ROTS.tsv", sep="\t", row.names=T, col.names=T, quote=F)






