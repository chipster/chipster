# ANALYSIS "Statistics"/"Sample size calculations with an adapted BH method" (Perform sample size calculations using an adapted Benjamini-Hochberg method.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv
# OUTPUT skewness.pdf, kurtosis.pdf, p-density.pdf, lambda.pdf, g.pdf, gamma.pdf, power.pdf, power.txt
# PARAMETER column METACOLUMN_SEL DEFAULT group (The phenodata column that divides the samples into exactly two groups.)
# PARAMETER assume.equal.variances [yes, no] DEFAULT no (Whether to treat the variances of the two groups as equal.)
# PARAMETER distribution [normal, student] DEFAULT normal (Whether to use the normal or the t distribution to calculate p values.)
# PARAMETER false.discovery.rate DECIMAL DEFAULT 0.1 (False discovery rate.)
# PARAMETER image.width INTEGER FROM 200 TO 6400 DEFAULT 2400 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 6400 DEFAULT 2400 (Height of the plotted network image)

# sample-size-with-bh.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2010-10-12

dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)
phenodata <- read.table("phenodata.tsv", header=TRUE, sep="\t")
groups <- phenodata[,column]

if(length(unique(groups[!is.na(groups)]))!=2)
   stop('CHIPSTER-NOTE: You need to have exactly two groups to run this analysis')

example.data = list(x=dat[,grep('chip', colnames(dat))], y=groups, genenames=rownames(dat), geneid=1:nrow(dat), samplelabels=phenodata$description)

# This program implements the method proposed in the article "Approximate sample size calculations
# with microarray data: an illustration", by J.A. Ferreira and A.H. Zwinderman, and references
# therein. The notation follows closely that of the article.
# This is a reworking of an earlier version written for S-plus. Among other things, this version
# uses t-statistics rather than z-statistics, so it requires the loading of the library QRMlib,
# which makes it possible to compute the modified Bessel function of the 3rd kind (in terms of
# which the characteristic function of the Student distribution is expressed: see Heyde and
# Leonenko (2005), Advances in Applied Probability, 37, 342-365).
# Date: 19/02/2008
#
################################################################################################################################################
# 1. Uploading and preliminary analysis of data
################################################################################################################################################
#library(PAMr)
#library(golubEsets)
#
#help(Golub_Merge)
#data(Golub_Merge)
library(moments)
#gene.expressions<-exprs(Golub_Merge)
#Golub.data<-list(x=matrix(gene.expressions,nrow=7129,ncol=72,dimnames=NULL),y=Golub_Merge$ALL.AML,
#genenames=row.names(gene.expressions),geneid=(1:7129),samplelabels=colnames(gene.expressions))
#
#example.data<-khan.data # khan.data, Golub.data
n<-nrow(example.data$x) 
Sample.of.test.statistics<-rep(0,n)
class.names<-unique(example.data$y)
samples1<-as.matrix(example.data$x[,example.data$y==class.names[1] & !is.na(example.data$y)]) # now works with samples
samples2<-as.matrix(example.data$x[,example.data$y==class.names[2] & !is.na(example.data$y)]) # with missing group information
N1<-ncol(samples1)
N2<-ncol(samples2)
N<-1/(1/N1+1/N2)
for(i in 1:n){
  prob <- TRUE
  try({
    Sample.of.test.statistics[i]<-t.test(samples1[i,],samples2[i,],alternative="two.sided",mu=0,paired=FALSE,var.equal=(assume.equal.variances=='yes'))$statistic
    prob <- FALSE
  }, silent=TRUE)
  if (prob)
    Sample.of.test.statistics[i] <- 0
}
#
#tests.of.normality<-matrix(0,nrow=n,ncol=2)
#sample.variances<-matrix(0,nrow=n,ncol=2)
sample.skewness<-matrix(0,nrow=n,ncol=2)
sample.kurtosis<-matrix(0,nrow=n,ncol=2)
for(i in 1:n){
  #sample.means<-c(mean(samples1[i,]),mean(samples2[i,]))
  #sample.variances[i,]<-c(var(samples1[i,]),var(samples2[i,]))
  #tests.of.normality[i,]<-c(shapiro.test(samples1[i,])$p.value,shapiro.test(samples2[i,])$p.value)
  sample.skewness[i,] <- c(skewness(samples1[i,]), skewness(samples2[i,]))
  sample.kurtosis[i,] <- c(kurtosis(samples1[i,]), kurtosis(samples2[i,]))
}
normal.skewness<-matrix(0,nrow=1000,ncol=2)
normal.kurtosis<-matrix(0,nrow=1000,ncol=2)
for (i in 1:1000) {
  rnd1 <- rnorm(ncol(samples1))
  rnd2 <- rnorm(ncol(samples2))
  normal.skewness[i,] <- c(skewness(rnd1), skewness(rnd2))
  normal.kurtosis[i,] <- c(kurtosis(rnd1), kurtosis(rnd2))
}
hist.norm <- function(x, norm = NULL, ...) {
  u <- seq(min(x, na.rm=TRUE)-1, max(x, na.rm=TRUE)+1, length=401)
  if (is.null(norm)) {
    d <- dnorm(u)
  } else {
    d <- dnorm(u, mean=mean(norm), sd=sd(norm))
  }
  h <- hist(x, plot=FALSE)
  hist(x, prob=TRUE, ylim=c(0,max(h$density, d)), ...)
  lines(u, d, col='red')
}
#par(mfrow=c(1,2))
#hist(tests.of.normality[,1],prob=TRUE,xlab="p-value of Shapiro test",main=as.character(class.names[1]))
#hist(tests.of.normality[,2],prob=TRUE,xlab="p-value of Shapiro test",main=as.character(class.names[2]))
#
#windows()
#aux<-c(min(log(sample.variances[,1]),log(sample.variances[,2])),max(log(sample.variances[,1]),log(sample.variances[,2])))
#plot(log(sample.variances[,1]),log(sample.variances[,2]),xlab=paste("log of variances in",as.character(class.names[1]),""),
#ylab=paste("log of variances in",as.character(class.names[2]),""),xlim=aux,ylim=aux)
#lines((aux[1]:aux[2]),(aux[1]:aux[2]),col=2)
#
#windows()
#aux1<-c(min(Sample.of.test.statistics),max(Sample.of.test.statistics))
#aux2<-c(min(log(sample.variances[,1]/sample.variances[,2])),max(log(sample.variances[,1]/sample.variances[,2])))
#plot(Sample.of.test.statistics,log(sample.variances[,1]/sample.variances[,2]),xlab="t-statistic",ylab="log of variance ratio",
#  xlim=aux1,ylim=aux2)
#
#for(j in 1:2){
#  aux1<-c(min(log(sample.variances[,1]/sample.variances[,2])),max(log(sample.variances[,1]/sample.variances[,2])))
#  aux2<-c(min(-log(tests.of.normality[,j])),max(-log(tests.of.normality[,j])))
#  plot(log(sample.variances[,1]/sample.variances[,2]),-log(tests.of.normality[,j]),xlab="log of variance ratio",
#    ylab="-log of p-value of Shapiro test",xlim=aux1,ylim=aux2,main=as.character(class.names[j]))
#  if(j==1){windows()}
#}
#
#for(j in 1:2){
#  aux1<-c(min(Sample.of.test.statistics),max(Sample.of.test.statistics))
#  aux2<-c(min(-log(tests.of.normality[,j])),max(-log(tests.of.normality[,j])))
#  plot(Sample.of.test.statistics,-log(tests.of.normality[,j]),xlab="t-statistic",ylab="-log of p-value of Shapiro test",xlim=aux1,ylim=aux2,
#    main=as.character(class.names[j]))
#  if(j==1){windows()}
#}
pdf(file='skewness.pdf', width=image.width/72, height=image.height/72)
par(mfrow=c(1,2))
hist.norm(sample.skewness[,1], normal.skewness[,1], main=as.character(class.names[1]), xlab='skewness')
hist.norm(sample.skewness[,2], normal.skewness[,2], main=as.character(class.names[2]), xlab='skewness')
dev.off()
pdf(file='kurtosis.pdf', width=image.width/72, height=image.height/72)
par(mfrow=c(1,2))
hist.norm(sample.kurtosis[,1], normal.kurtosis[,1], main=as.character(class.names[1]), xlab='kurtosis')
hist.norm(sample.kurtosis[,2], normal.kurtosis[,2], main=as.character(class.names[2]), xlab='kurtosis')
dev.off()
# Reducing the data by considering only those genes with "good" goodness-of-fit results (those with p-values of the Shapiro test
# above a certain small value)
#n
#small.value<-70/n
#example.data<-reduce.data(small.value,example.data,tests.of.normality)
#nrow(example.data$x)
#
# Auxiliary function:
#reduce.data<-function(small.value,data,tests.of.normality)
#{
#  new.p.values<-1-(1-apply(tests.of.normality,1,min))^2
#  ratios.of.variances<-sample.variances[,1]/sample.variances[,2]
#  aux1<-(1:length(new.p.values))
#  #indices<-aux1[new.p.values>=small.value]
#  indices<-aux1[ratios.of.variances>=2/3 & ratios.of.variances<=3/2 & new.p.values>=small.value]
#  reduced.data<-list(x=data$x[indices,],y=data$y,genenames=data$genenames[indices],
#    geneid=data$geneid[indices],samplelabels=data$samplelabels,batchlabels=data$batchlabels)
#  return(reduced.data)
#}
#
################################################################################################################################################
# 2. Computing the sample of p-values, and plotting a density estimate of the p-values
################################################################################################################################################
# Auxiliary function:
beta.kernel.density.at.x<-function(x,data,a.n){
  return(ifelse(min(x,1-x)<=0,0,mean(dbeta(data,(x/a.n)+1,((1-x)/a.n)+1))))
}
#
if (distribution == 'student') {
    Sample.of.p.values<-2*(1-pt(abs(Sample.of.test.statistics),N1+N2-2,0)) # student
} else {
    Sample.of.p.values<-2*(1-pnorm(abs(Sample.of.test.statistics),0)) # normal
}
#
resolution<-100
c<-1
a.n<-n^((c-1)/2) #Smoothing parameter 
x<-seq(1/resolution,1-1/resolution,1/resolution)
density.estimate.of.p.values<-sapply(x,beta.kernel.density.at.x,Sample.of.p.values,a.n)
pdf(file='p-density.pdf', width=image.width/72, height=image.height/72)
plot(x,density.estimate.of.p.values,type='l',xlab="p-value",ylab="density",col=1,xlim=c(0,1),
ylim=c(0,1.1*max(density.estimate.of.p.values)))
rough.estimate.of.gamma<-beta.kernel.density.at.x(1-0.0001,Sample.of.p.values,a.n)
rough.estimate.of.gamma
lines(x,rep(rough.estimate.of.gamma,length(x)),col=2)
dev.off()
#
################################################################################################################################################
# 3. Definition of the integrands involved in the calculation of the estimate of lambda
################################################################################################################################################
# Auxiliary function:
Student.characteristic.function<-function(t,nu){
   aux<-ifelse(t==0,0,log(besselK(abs(t)*sqrt(nu),nu/2,expon.scaled=FALSE)*((abs(t)*sqrt(nu))^(nu/2))*(2^(1-nu/2)))-lgamma(nu/2))
   return(exp(aux))
}
#
first.integrand.at.t<-function(t,theta,statistics,n,a.n,N,nu){   
   cosine.t.theta<-cos(t*theta*sqrt(N))
   sine.t.theta<-sin(t*theta*sqrt(N))
   if (distribution == 'student') {
      Phi.star.t<-max(Student.characteristic.function(t,nu),0.000000000000001) # student, can sometimes be zero
   } else {
      Phi.star.t<-exp(-(t^2)/2) # normal
   }
#  Note the two choices available for Kappa.star.a.n.t
#  Kappa.star.a.n.t<-max(0,1.0-abs(a.n*t))
   Kappa.star.a.n.t<-pmax(0,1.0-abs(a.n*t)^2)^3
   Re.M.star.t<-rep(0,length(t))
   Im.M.star.t<-rep(0,length(t))
   for(i in 1:length(t)){
   Re.M.star.t[i]<-mean(cos(t[i]*statistics))
   Im.M.star.t[i]<-mean(sin(t[i]*statistics))
   }
   first.integrand<-Kappa.star.a.n.t*(cosine.t.theta*Re.M.star.t+sine.t.theta*Im.M.star.t)/Phi.star.t
   return(first.integrand)
}
#
second.integrand.at.t<-function(t,theta,statistics,n,a.n,N){   
   cosine.t.theta<-cos(t*theta*sqrt(N))
#  Note the two choices available for Kappa.star.a.n.t
#  Kappa.star.a.n.t<-max(0,1.0-abs(a.n*t))
   Kappa.star.a.n.t<-pmax(0,1.0-abs(a.n*t)^2)^3
   second.integrand<-Kappa.star.a.n.t*cosine.t.theta
   return(second.integrand)
}
#
# Examples of calls:
nu<-N1+N2-2
M<-1000
d.n<-1/(a.n*M)
t<-seq(from=-d.n*M,by=d.n,to=d.n*M)
a.n<-1.0/sqrt(log(n)) #a.n is the smoothing parameter used to smooth the empirical c.f.
first.integrand.at.t(t,0,Sample.of.test.statistics,n,a.n,N,nu)
second.integrand.at.t(t,0,Sample.of.test.statistics,n,a.n,N)
#
################################################################################################################################################
# 4. Computing and plotting estimates of lambda; each estimate is obtained by choosing a value of "initial.gamma"
################################################################################################################################################
# Auxiliary function:
estimate.of.lambda.at.theta<-function(theta,a.n,statistics,n,N,nu){
  first.integral<-integrate(first.integrand.at.t,lower=-1/a.n,upper=1/a.n,theta=theta,
  subdivisions=100000,abs.tol=0.00001,statistics=statistics,n=n,a.n=a.n,N=N,nu=nu)$value
  second.integral<-integrate(second.integrand.at.t,lower=-1/a.n,upper=1/a.n,theta=theta,
  subdivisions=100000,abs.tol=0.00001,statistics=statistics,n=n,a.n=a.n,N=N)$value
  return(c(first.integral,second.integral))
}
#
nu<-N1+N2-2
step.theta<-0.01
theta.vector<-seq(-6,6,step.theta) #Points at which lambda is to be estimated
a.n<-0.5*sqrt(nu)/log(n) #a.n is the smoothing parameter used to smooth the empirical c.f.
statistics<-Sample.of.test.statistics
initial.gamma<-0.7
#
# NB: "components.lambda" needs to be computed only ONCE; it is to be used in everything that follows
components.lambda<-t(apply(as.matrix(theta.vector),1,estimate.of.lambda.at.theta,a.n,statistics,n,N,nu))
estimate.lambda.of.theta<-sqrt(N)*(components.lambda[,1]-initial.gamma*components.lambda[,2])/((1-initial.gamma)*(2*pi))
estimate.lambda.of.theta<-pmax(rep(0,length(theta.vector)),estimate.lambda.of.theta)
Const<-sum(estimate.lambda.of.theta)*step.theta
estimate.lambda.of.theta<-estimate.lambda.of.theta/Const
#Plotting estimate.lambda.of.theta
pdf(file='lambda.pdf', width=image.width/72, height=image.height/72)
plot(theta.vector,estimate.lambda.of.theta,type='l',col=3,ylim=c(0,1.2*max(estimate.lambda.of.theta)),xlab=quote(theta),
ylab=quote(lambda(theta)))
title("Estimated density of effect sizes")
dev.off()
#
################################################################################################################################################
# 5. Computing and comparing G.n and G.hat.n
################################################################################################################################################
# First auxiliary function:
edf.at.a.point<-function(x,random.sample)
{
  sample.size<-length(random.sample)
  edf.at.x<-length(random.sample[random.sample<=x])/sample.size
  return(edf.at.x)
}
#
# Second auxiliary function:
Big.gamma.at.theta<-function(theta,u,nu,N)
{
  if (distribution == 'student') {
    return(1-pt(qt(1-u/2,nu)-theta*sqrt(N),nu)+pt(-qt(1-u/2,nu)-theta*sqrt(N),nu)) # student
  } else {
    return(1-pnorm(qnorm(1-u/2)-theta*sqrt(N))+pnorm(-qnorm(1-u/2)-theta*sqrt(N))) # normal
  }
}
#
# Third auxiliary function:
G.hat.n.at.a.point<-function(u,theta.vector,step.theta,estimate.lambda.of.theta,nu,N)
{
  aux<-sapply(as.matrix(theta.vector),Big.gamma.at.theta,u,nu,N)
  return(aux%*%estimate.lambda.of.theta*step.theta)
}
#
initial.gamma<-0.54
step.u<-0.01
u<-seq(0,1,step.u)
#
estimate.lambda.of.theta<-sqrt(N)*(components.lambda[,1]-initial.gamma*components.lambda[,2])/((1-initial.gamma)*(2*pi))
estimate.lambda.of.theta<-pmax(rep(0,length(theta.vector)),estimate.lambda.of.theta)
Const<-sum(estimate.lambda.of.theta)*step.theta
estimate.lambda.of.theta<-estimate.lambda.of.theta/Const
#
G.n.hat.u<-sapply(as.matrix(u),G.hat.n.at.a.point,theta.vector,step.theta,estimate.lambda.of.theta,nu,N)
H.n.u<-sapply(u,edf.at.a.point,Sample.of.p.values)
G.n.tilde.u<-(H.n.u-initial.gamma*u)/(1-initial.gamma)
pdf(file='g.pdf', width=image.width/72, height=image.height/72)
plot(u,G.n.tilde.u,type='l',xlab=quote(u),ylab=quote(G(u)),lty=c(1,3),col=c(3,2),xlim=c(0,1),ylim=c(0,max(G.n.tilde.u)))
lines(u,G.n.hat.u,col=2)
legend(c(0.6,0.8),c(0.6,0.8),c("non-parametric","semi-parametric"),lty=c(1,3),col=c(3,2),bty="n")
title("Estimates of G")
dev.off()
distance<-sum(abs(G.n.hat.u-G.n.tilde.u))*step.u
#
################################################################################################################################################
# 6. Determination of the value of Gamma.hat that minimizes D.n
################################################################################################################################################
# Auxiliary function:
distance.between.estimates.of.G.n<-function(gamma.hat,u,step.u,theta.vector,step.theta,components.lambda.of.theta,nu,N)
{
  estimate.lambda.of.theta<-sqrt(N)*(components.lambda[,1]-gamma.hat*components.lambda[,2])/((1-gamma.hat)*(2*pi))
  estimate.lambda.of.theta<-pmax(rep(0,length(theta.vector)),estimate.lambda.of.theta)
  Const<-sum(estimate.lambda.of.theta)*step.theta
  estimate.lambda.of.theta<-estimate.lambda.of.theta/Const
  G.n.hat.u<-sapply(as.matrix(u),G.hat.n.at.a.point,theta.vector,step.theta,estimate.lambda.of.theta,nu,N)
  H.n.u<-sapply(u,edf.at.a.point,Sample.of.p.values)
  G.n.tilde.u<-(H.n.u-gamma.hat*u)/(1-gamma.hat)
  distance<-sum(abs(G.n.hat.u-G.n.tilde.u))*step.u
  return(distance)
}
#
step.gamma<-0.01
gamma.vector<-seq(step.gamma,1-step.gamma,step.gamma)
distances<-sapply(as.matrix(gamma.vector),distance.between.estimates.of.G.n,u,step.u,theta.vector,step.theta,components.lambda.of.theta,nu,N)
#
minimum.distance<-min(distances)
gamma.hat<-gamma.vector[distances<=minimum.distance]
#
pdf(file='gamma.pdf', width=image.width/72, height=image.height/72)
plot(gamma.vector,distances,type='l',xlab=quote(~~gamma),ylab=quote(D[n](gamma)),col=1,ylim=c(0,1))
dev.off()
#
################################################################################################################################################
# 7. Computation of the sample size required to achieve a certain power subject to a given false discovery rate (f.d.r.)
################################################################################################################################################
Gamma<-gamma.hat
#delta<-0.1 # nominal f.d.r.
delta<-false.discovery.rate
r<-delta*(1-Gamma)/(Gamma*(1-delta))
#
step.u<-0.01
u<-seq(0,1,step.u)
#Nprime<-N
#
estimate.lambda.of.theta<-sqrt(N)*(components.lambda[,1]-initial.gamma*components.lambda[,2])/((1-initial.gamma)*(2*pi))
estimate.lambda.of.theta<-pmax(rep(0,length(theta.vector)),estimate.lambda.of.theta)
Const<-sum(estimate.lambda.of.theta)*step.theta
estimate.lambda.of.theta<-estimate.lambda.of.theta/Const
#
samplesize.coefs <- c(0.25, 0.5, 0.75, 1:30)
samplesize.N <- 1/(1/N1/samplesize.coefs + 1/N2/samplesize.coefs)
samplesize <- samplesize.coefs * (N1+N2)
if (samplesize[length(samplesize)-1] > 300) {
  last <- max(which(samplesize <= 300)) + 1
  samplesize.N <- samplesize.N[1:last]
  samplesize <- samplesize[1:last]
}
power <- numeric(length=length(samplesize.N))
names(power) <- samplesize.N
for(Nprime in samplesize.N) {
    G.n.hat.u<-sapply(as.matrix(u),G.hat.n.at.a.point,theta.vector,step.theta,estimate.lambda.of.theta,nu,Nprime)
    #
    u.star<-max(u[G.n.hat.u>u/r])
    this.power<-u.star/r # laske tama
    power[as.character(Nprime)] <- this.power
    #plot(u,G.n.hat.u,type='l',xlab=quote(u),ylab=quote(G(u)),lty=c(1,3),col=c(3,2),xlim=c(0,u.star*1.5),ylim=c(0,max(G.n.tilde.u)))
    #lines(u,u/r)
    #title(paste("Power calculation for N'=",round(Nprime,digits=1),""))
}
pdf(file='power.pdf', width=image.width/72, height=image.height/72)
plot(samplesize, power, type='b', main=paste('Power calculation for FDR = ', round(100*delta), '%', sep=''),
  xlim=c(0, 300), ylim=0:1, ylab='Average power',
  xlab=paste('Number of samples (', round(100*N1/(N1+N2)), '% ', class.names[1], ', ', round(100*N2/(N1+N2)), '% ', class.names[2], ')', sep=''))
dev.off()
#
write.table(data.frame(samplesize, power), file="power.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')
#
################################################################################################################################################

# EOF
