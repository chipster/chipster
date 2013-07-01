# TOOL plot-dendrogram.R: Dendrogram (Creates a dendrogram of samples using normalized data with Pearson correlation and average linkage method. The branches of the tree are colored according to the selected number of groups. Clustering is done using the function hcluster in which the parameter correlation envokes computation of pearson type of distances.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT dendrogram-color.pdf: dendrogram-color.pdf 
# OUTPUT dendrogram-bw.pdf: dendrogram-bw.pdf 
# PARAMETER cluster: cluster TYPE [genes: genes, chips: chips] DEFAULT chips (What to cluster)
# PARAMETER column: column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to color next to the tree)
# PARAMETER number.of.groups: number.of.groups TYPE INTEGER FROM 2 TO 20 DEFAULT 2 (How many groups to color to the tree)
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Dendrogram
# JTT 3.10.2007
#
# MG 25.11.2010
# Increased the gene/sample limit to 20000

# Renaming variables
w<-image.width
h<-image.height
gr<-number.of.groups 
margin<-cluster

# Loading the libraries
library(fpc)
library(A2R)
library(amap)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Loads the phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
groups<-phenodata[,pmatch(column,colnames(phenodata))]

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]
colnames(dat2)<-gsub(" ", "", phenodata$description)

# Manipulate data depending on what to cluster
if(margin=="chips") {
dat2<-t(dat2) 
}

if (nrow(dat2) > 20000) {
stop("Hierarchical clustering dendrogram can be plotted on maximum 2000 of genes/samples");
}

# Does the clustering
clust<-hcluster(x=dat2, method="correlation", link="average")
ct<-cutree(clust, gr)


# A2Rplot code from http://addictedtor.free.fr/packages/A2R/lastVersion/R/A2R
.packageName <- "A2R"
arimaSelect <- function(
  serie,                # the time serie to fit
  order,                # see stats::arima
  include.mean = FALSE, # idem
  alpha = 0.05          # risk
  ){

   n_max <- sum( order[ c(1,3) ] )
   n     <- length(serie)
   
   #! the matrix of coeffs, each row for each step
   coeff     <- matrix( NA , nrow = n_max , ncol = n_max )
   
   #! same with t-stat p-value's
   pval      <- coeff
   aic       <- rep(NA, n_max)
   liste     <- rep( list(NULL) , n_max )
   fixed     <- aic
   
   go <- TRUE
   i  <- 1
   while(go){
     
     arima.out <- arima(serie,
                        order        = order,
                        include.mean = include.mean,
                        method       = "ML", 
                        fixed        = fixed)
     liste[[i]] <- arima.out
     
     .coeff <- arima.out$coef
     .sd    <- rep(0, n_max)
     .sd[arima.out$mask] <- sqrt(diag(arima.out$var.coef))

     .pval  <- round( 2 * (1 - pt(abs( .coeff / .sd  ),
                                  df = n - n_max - 1 ) ) , 5)

     coeff[i,] <- .coeff
     pval[i,]  <- .pval
     aic[i]    <- arima.out$aic
     fixed[which.max(.pval)] <- 0
     go     <- any(.pval > alpha, na.rm = TRUE)
     i      <- i+1

     
   }
   
   
res <- list(coeff     = coeff,
            pval      = pval,
            aic       = aic, 
            listArima = liste)
class(res) <- "arimaSelect"
res
            
}

plot.arimaSelect <- function(x, choix, ...){
	noms <- names(x$listArima[[1]]$coef)
	coeff <- x$coeff
	k <- min(which(is.na(coeff[,1])))-1
	coeff <- coeff[1:k,]
	pval  <- x$pval[1:k,]
	aic   <- x$aic[1:k]
	coeff[coeff==0] <- NA
	n <- ncol(coeff)
	if(missing(choix)) choix <- k
	
	layout(matrix(c(1,1,1,2,
	                3,3,3,2,
			3,3,3,4,
			5,6,7,7),nr=4),
               widths=c(10,35,45,15),
	       heights=c(30,30,15,15))
	couleurs <- rainbow(75)[1:50]#(50)
	ticks <- pretty(coeff)
	
	### graph AIC
	par(mar=c(1,1,3,1))
	plot(aic,k:1-.5,type="o",pch=21,bg="blue",cex=2,axes=FALSE,lty=2,xpd=NA)
	points(aic[choix],k-choix+.5,pch=21,cex=4,bg=2,xpd=NA)
	#axis(3)
	title("aic",line=2)
	
	par(mar=c(3,0,0,0))	
	plot(0,axes=FALSE,xlab="",ylab="",xlim=range(ticks),ylim=c(.1,1))
	rect(xleft  = min(ticks) + (0:49)/50*(max(ticks)-min(ticks)),
	     xright = min(ticks) + (1:50)/50*(max(ticks)-min(ticks)),
	     ytop   = rep(1,50),
	     ybottom= rep(0,50),col=couleurs,border=NA)
	axis(1,ticks)
	rect(xleft=min(ticks),xright=max(ticks),ytop=1,ybottom=0)
	text(mean(coeff,na.rm=T),.5,"coefficients",cex=2,font=2)
	
	
	par(mar=c(1,1,3,1))
	image(1:n,1:k,t(coeff[k:1,]),axes=FALSE,col=couleurs,zlim=range(ticks))
	for(i in 1:n) for(j in 1:k) if(!is.na(coeff[j,i])) {
		if(pval[j,i]<.01)                            symb = "green"
		else if( (pval[j,i]<.05) & (pval[j,i]>=.01)) symb = "orange"
		else if( (pval[j,i]<.1)  & (pval[j,i]>=.05)) symb = "red"
		else                                         symb = "black"
		polygon(c(i+.5   ,i+.2   ,i+.5   ,i+.5),
		        c(k-j+0.5,k-j+0.5,k-j+0.8,k-j+0.5),
			col=symb)
		
		#points(i+.4,k-j+.6,pch=21,bg=symb)
		#text(i+.5,k-j+.8,round(pval[j,i],2),pos=2,cex=.8)
		if(j==choix)  {
			rect(xleft=i-.5,
			     xright=i+.5,
			     ybottom=k-j+1.5,
			     ytop=k-j+.5,
			     lwd=4)
			text(i,
			     k-j+1,
			     round(coeff[j,i],2),
			     cex=1.2,
			     font=2)
		}
		else{
			rect(xleft=i-.5,xright=i+.5,ybottom=k-j+1.5,ytop=k-j+.5)
			text(i,k-j+1,round(coeff[j,i],2),cex=1.2,font=1)
		}
	}
	axis(3,1:n,noms)


	par(mar=c(0.5,0,0,0.5))	
	plot(0,axes=FALSE,xlab="",ylab="",type="n",xlim=c(0,8),ylim=c(-.2,.8))
	cols <- c("green","orange","red","black")
	niv  <- c("0","0.01","0.05","0.1")
	for(i in 0:3){
		polygon(c(1+2*i   ,1+2*i   ,1+2*i-.5   ,1+2*i),
		        c(.4      ,.7      , .4        , .4),
			col=cols[i+1])
		text(2*i,0.5,niv[i+1],cex=1.5)	
		}
	text(8,.5,1,cex=1.5)
	text(4,0,"p-value",cex=2)
	box()
	
	residus <- x$listArima[[choix]]$res
	
	par(mar=c(1,2,4,1))
	acf(residus,main="")
	title("acf",line=.5)
	
	par(mar=c(1,2,4,1))
	pacf(residus,main="")
	title("pacf",line=.5)
	
	par(mar=c(2,2,4,1))
	qqnorm(residus,main="")
	title("qq-norm",line=.5)
	
}

"._a2r_hclu"       <- NULL # to receive an hclust object when 
                           # A2Rplot.hclust is called

"._a2r_counter"       <- NA # a counter used in A2Rplot.hclust
"._a2r_height_cut"    <- NA

"._a2r_envir"         <- NA
"._a2r_group"         <- NA


#===============================================================================
"A2Rplot" <- function(x,...){
  UseMethod("A2Rplot")
}
#===============================================================================
"A2Rplot.default" <- function(x,...){
  plot(x,...)
}
#===============================================================================
"A2Rplot.hclust" <- function(
  x ,             # an hclust object to draw
  k        = 2,   # the number of groups
  col.up   = "black",
  col.down = rainbow(k),
  lty.up   = 2,
  lty.down = 1,
  lwd.up   = 1,
  lwd.down = 2,
  type     = c("rectangle","triangle"),
  knot.pos = c("mean","bary","left","right","random"),
  criteria,
  fact.sup,
  show.labels=TRUE,
  only.tree=FALSE,
  main     = paste("Colored Dendrogram (",k," groups)"),
  boxes    = TRUE,
  members,
  ...
){

  if(missing(members)) members <- NULL
  opar <- par(no.readonly=TRUE)
  knot.pos <- match.arg(knot.pos)
  type     <- match.arg(type)
  # tests
  if(k<2) 
    stop("k must be at least 2")  
    
  ._a2r_counter    <<- 0
  ._a2r_hclu       <<- x

  ._a2r_envir      <<- environment()
  nn <- length(x$order) - 1

  ._a2r_height_cut <<- mean(x$height[nn-k+1:2])
  ._a2r_group      <<- 0
  
  n.indiv   <- length(x$order)
  groups.o  <- cutree.order(x, k=k)[x$order]
  
  bottom <- if(is.null(members)) 0 else x$height[nn] * -.2 
  
  if(only.tree){
    if(is.null(members)) plot(0,type="n",xlim=c(0.5,n.indiv+.5), ylim=c(bottom,x$height[nn]), xaxs="i", axes=FALSE, xlab="",ylab="") 
    else                 plot(0,type="n",xlim=c(0.5,sum(members)+.5), ylim=c(bottom,x$height[nn]), xaxs="i", axes=FALSE, xlab="",ylab="")
    #call to the ** recursive function ** .rec.hclust
    .rec.hclust(nn, col=col.up, lty=lty.up, lwd=lwd.up)
    
    if(boxes){
      axis(2)
      box()
    }
    return(NULL)
  }
  
  # prepare the layout
  matlayout <- matrix(c(2,4,6,1,3,5), nc=2, nr=3)
  widths    <- c(1,9)
  heights   <- c(8,1,1)
  if(!show.labels){
      matlayout <- matrix(c(2,4,1,3), nc=2, nr=2)
      widths    <- c(1,9)
      heights   <- c(9,1)
  }
  if(!missing(fact.sup) ) {
    heights   <- c(8,1,1)
  }
  if(missing(criteria) & missing(fact.sup)){
    matlayout <- matrix(c(2,4,1,3), nc=2, nr=2)
      widths    <- c(1,9)
      heights   <- c(9,1)
    
  }
  layout(matlayout, width=widths, height=heights)
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ The tree (1)
  par(mar=c(0,0,3,4))
  if(is.null(members)) plot(0,type="n",xlim=c(0.5,n.indiv+.5), ylim=c(bottom,x$height[nn]), xaxs="i", axes=FALSE, xlab="",ylab="") 
  else plot(0,type="n",xlim=c(0.5,sum(members)+.5), ylim=c(bottom,x$height[nn]), xaxs="i", axes=FALSE, xlab="",ylab="") 
  #call to the ** recursive function ** .rec.hclust
  .rec.hclust(nn, col=col.up, lty=lty.up, lwd=lwd.up)
  title(main)
  if(boxes){
    box()
    axis(4)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Criteria (2)
  if(!missing(criteria)){
    par(mar=c(0,0,3,0))
    plot(0,
         type="n",
         xlim=range(criteria), 
         ylim=c(0,x$height[nn]), 
         axes=FALSE, 
         xlab="",
         ylab="")
    par(las=2)
    n.crit <- length(criteria)
    heights.cut <- ( tail(x$height,n.crit) + 
                     tail(x$height,n.crit+1)[-(n.crit+1)] ) / 2
    heights.cut <- rev(heights.cut)
                   
    points(criteria   , heights.cut   , pch=21, bg="red", type="o")
    points(criteria[k-1], heights.cut[k-1], pch=21, cex=2, bg="blue", xpd=NA)
    if(boxes){
      axis(3)
      box()
    }
  }
  else{
    par(mar=c(0,0,3,0))
    plot(0,
         type="n",
         xlim=c(0,1), 
         ylim=c(0,1), 
         axes=FALSE, 
         xlab="",
         ylab="")
  }

  if(show.labels){
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Name of the observations (3)
    par(mar=c(0,0,0,4))
    par(srt=90)
    obs.labels <- toupper(substr(x$labels[x$order],1,6))
    if(is.null(members)) {
      plot(0,type="n",xlim=c(0.5,n.indiv+.5), ylim=c(0,1), xaxs="i", axes=FALSE, xlab="",ylab="") 
      text(1:n.indiv            , 0, obs.labels, pos=4, col=col.down[groups.o])
    }
    else{
      plot(0,type="n",xlim=c(0.5,sum(members)+.5), ylim=c(0,1), xaxs="i", axes=FALSE, xlab="",ylab="") 
      xo <-   members[x$order]
      text(cumsum(xo)-xo/2, 0, obs.labels, pos=4, col=col.down[groups.o])
    }
    par(srt=0)
    if(boxes){
      box()
    }
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Labels (4)
    par(mar=c(0,0,0,0))
    plot(0,type="n",xlim=c(0,1), ylim=c(0,1), xaxs="i", axes=FALSE, xlab="",ylab="") 
    text(.5,.5,"Labels")
    if(boxes){
      box()
    }
      
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Quali (5,6)
  if(!missing(fact.sup)){
    quali  <- as.factor(fact.sup)[x$order]
    quanti <- as.numeric(quali)

    par(mar=c(1,0,0,4))
    n.levels <- length(levels(quali))
    plot(0,type="n",
         xlim=c(0.5,n.indiv+.5), 
         ylim=c(0,n.levels), 
         xaxs="i", yaxs="i",axes=FALSE, xlab="",ylab="") 
        
    rect(xleft    = (1:n.indiv)-.5,
         xright   = (1:n.indiv)+.5,
         ybottom  = quanti-1, 
         ytop     = quanti,
         col      = col.down[groups.o])
    par(las=1)
    axis(4, (1:n.levels)-.5,levels(quali), tick=FALSE)
      
    if(boxes){
      box()
    }
    
    
    par(mar=c(1,0,0,0))
    plot(0,type="n",xlim=c(0,1), ylim=c(0,1), xaxs="i", axes=FALSE, xlab="",ylab="") 
    text(.5,.5,deparse(substitute(fact.sup)))
    if(boxes){
      box()
    }
  }
  
  
  par(opar) # reset parameter
}

#===============================================================================

".rec.hclust" <- function(
  index, # index of the current tree to draw
  lwd = 1,
  lty = 1,
  col = "black"){

  members <- get('members', envir= ._a2r_envir) 
  bottom  <- get('bottom',  envir= ._a2r_envir) 
  if(index<0){ # it is a leaf
    if(is.null(members)){
       ._a2r_counter <<- ._a2r_counter + 1
       return(list( x = ._a2r_counter,
                    n = 1))       
    }
    else{
      cc <- ._a2r_counter
      mm <- members[-index]
      polygon(x  = c(cc, cc+mm/2, cc+mm),
              y  = c(bottom, 0, bottom),
              col= col, 
              border = col, 
              lwd=lwd)
      ._a2r_counter <<- ._a2r_counter + mm
      return(list(x = cc+mm/2,
                  n = mm))
    }
  }
  
  h.m   <- ._a2r_hclu$height[index]

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~ do left
  index.l  <- ._a2r_hclu$merge[index,1]
  
  h.l <- if(index.l<0) 0 else ._a2r_hclu$height[index.l]
  if(h.l<._a2r_height_cut & h.m > ._a2r_height_cut){
      ._a2r_group <<- ._a2r_group + 1
      col.l <- get("col.down",envir=._a2r_envir)[._a2r_group]
      lwd.l <- get("lwd.down",envir=._a2r_envir)
      lty.l <- get("lty.down",envir=._a2r_envir)
  }
  else{
      col.l <- col
      lwd.l <- lwd
      lty.l <- lty
  }
  out.l   <- .rec.hclust(index.l, col=col.l, lty=lty.l, lwd=lwd.l)
  x.l     <- out.l$x
  n.l     <- out.l$n
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ do right
  index.r  <- ._a2r_hclu$merge[index,2]
  h.r <- if(index.r<0) 0 else ._a2r_hclu$height[index.r]
  if(h.r<._a2r_height_cut & h.m > ._a2r_height_cut){
      ._a2r_group <<- ._a2r_group + 1
      col.r <- get("col.down",envir=._a2r_envir)[._a2r_group]
      lwd.r <- get("lwd.down",envir=._a2r_envir)
      lty.r <- get("lty.down",envir=._a2r_envir)
  }
  else{
      col.r <- col
      lwd.r <- lwd
      lty.r <- lty
  }
  out.r   <- .rec.hclust(index.r, col=col.r, lty=lty.r, lwd=lwd.r)
  x.r     <- out.r$x
  n.r     <- out.r$n
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ draw what you have to draw
  
  type <- get("type",envir=._a2r_envir)
  x.m  <- (x.r + x.l) / 2  
  n    <- n.r + n.l
  x.b  <- (n.r * x.r + n.l * x.l) / n

  
  knot.pos <- get("knot.pos",envir=._a2r_envir) 
  
  x <- switch(knot.pos,
          mean = x.m,
          left = x.l,
          right= x.r,
          random = x.l + runif(1)*(x.r-x.l),
          bary   = x.b)

          
          
  if(type=="rectangle"){
    segments(x0  = c(x.l, x.l, x.r),
             x1  = c(x.l, x.r, x.r),
             y0  = c(h.l, h.m, h.r),
             y1  = c(h.m, h.m, h.m),
             col = col,
             lty = lty,
             lwd = lwd)
  }
  if(type =="triangle"){
    segments(x0  = c(x.l, x.r),
             x1  = c(x  , x),
             y0  = c(h.l, h.r),
             y1  = c(h.m, h.m),
             col = col,
             lty = lty,
             lwd = lwd)
  }
          
          
  list(x=x,n=n)
}
#===============================================================================
"cutree.order" <- function(hclu, k=NULL, h=NULL){
  
  coupe <- cutree(hclu,k=k, h=h)

  coupe.or <- coupe[hclu$order]
  coupe.out<- rep(NA,length(coupe))
  j <-  1 #
  k <-  coupe.or[1]
  for(i in 1:length(coupe)){
    if(coupe.or[i]==k) next
    else{
      coupe.out[which(coupe==k)] <- j
      j <- j + 1
      k <- coupe.or[i]
    }
  }
  coupe.out[is.na(coupe.out)] <- j
  names(coupe.out) <- names(coupe)
  coupe.out
}
#===============================================================================

"spie" <- function(firstPartition, secondPartition){
  
  if(length(firstPartition)!=length(secondPartition)) 
    stop("'firstPartition' and 'secondPartition' have different lengths")
  
  if(sum(firstPartition)!=sum(secondPartition))
    stop("'firstPartition' and 'secondPartition' doesn't sum the same")
  
  angles <- cumsum(c(0, 2 * pi * firstPartition / sum(firstPartition)))
  
  radii  <- sqrt( (secondPartition / sum(secondPartition)) /
                  (firstPartition  / sum(firstPartition)) )
  
  namesSlices <- if(is.null(names(firstPartition))) 
                   1:length(firstPartition) 
                 else 
                   names(firstPartition)  
  
  structure(list(angles=angles, 
                 radii=radii, 
                 firstPartition=firstPartition,
                 secondPartition=secondPartition, 
                 namesSlices = namesSlices), 
            class="spie")
}

"plot.spie" <- 
  function(
  x, 
  multi, 
  col = rainbow(length(x$radii)),
  ...){
  
    maxRadii <- max(x$radii)
  
  grid.newpage()
  
  pushViewport(viewport(layout=grid.layout(1,1,respect=TRUE)))
  pushViewport(dataViewport(maxRadii*c(-1.1,1.1), 
                            maxRadii*c(-1.1,1.1),
                            layout.pos.col=1,
                            layout.pos.row=1))
    
  
  if(!missing(multi)){
    
      grid.circle(x=0, y=0, r=sqrt(multi), gp=gpar(col="gray"), default.units="native")
    
  
  }  
  
  
  for(i in 1:length(x$radii)){
    
    theta <- seq(x$angles[i], x$angles[i+1], length=100)
    
    grid.polygon(x   = c(0,  cos(theta) ,0),
                 y   = c(0,  sin(theta) ,0) , 
                 gp  = gpar(fill=col[i]),
                 default.units="native")    

    grid.polygon(x   = c(0, x$radii[i] * cos(theta) ,0),
                 y   = c(0, x$radii[i] * sin(theta) ,0) , 
                 gp  = gpar(fill=col[i], lwd=2),
                 default.units="native")    
           
                 
    angleAnn <- mean(x$angles[i+0:1])
    maxx <- max(1, x$radii[i])+maxRadii/10
    
    
    grid.rect( x = cos(angleAnn)*maxx,
               y = sin(angleAnn)*maxx, 
               width = 1.5*stringWidth(x$namesSlices[i]),
               height = 1.5*stringHeight(x$namesSlices[i]),
               default.units="native",
               gp = gpar(col=col[i], fill="white", lwd=2))

    
    grid.text(x$namesSlices[i], 
              x=cos(angleAnn)*maxx, 
              y=sin(angleAnn)*maxx,
              default.units="native")
    
  }
  
  
  
  
  if(!missing(multi)){ 
    
      grid.lines(x=unit(0,"native"), 
                 y=unit(c(0, max(sqrt(multi))), "native"), gp=gpar(col="gray"))
    
    for(i in multi){
      st <- paste("x", i)
      sw <- stringWidth(st)
      sh <- stringHeight(st)
      

      grid.rect( x = unit(0, "native"),
               y = unit(sqrt(i),"native"), 
               width = 1.5*sw,
               height = 1.5*sh , gp=gpar(fill="white", col="gray"))
               
      grid.text( st , 0, sqrt(i), default.units="native")
      
    }
    
  }
  
  upViewport(2)
 
}


# Plotting
pdf(file="dendrogram-color.pdf", width=w/72, height=h/72)
if(margin=="chips") {
A2Rplot(clust, k=gr, fact.sup=groups) 
} else {
A2Rplot(clust, k=gr) 
}
dev.off()

pdf(file="dendrogram-bw.pdf", width=w/72, height=h/72)
plot(clust, hang=0.1)
dev.off()

