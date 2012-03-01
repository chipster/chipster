library(CGHcall)
library(CGHregions)

setMethod("plot", signature(x="cghRaw", y="missing"),
function (x, y, dotres=1, ylimit=c(-2,5), ylab=expression(log[2]~ratio), build="GRCh37",... )
{
    for (i in 1:ncol(x)) {
        cat("Plotting sample", sampleNames(x)[i], "\n")
        chrom           <- chromosomes(x)
        data            <- data.frame(chrom, bpstart(x), copynumber(x)[,i])
        colnames(data)  <- c("chromosome", "position", "ratio")
        chrom.labels    <- as.character(unique(chrom))
        pos             <- bpstart(x)
        chrom.ends <- .getCumulativeChromosomeEnds(build)[1:max(chrom)]
        for (j in 2:max(chrom))
            pos[chrom == j] <- pos[chrom == j] + chrom.ends[j-1]
        nclone <- length(chrom)
        whichtoplot <- seq(1,nclone,by=dotres) #added 15/06/2009
        plot(pos[whichtoplot], data[whichtoplot,3], cex=.1, main=sampleNames(x)[i], ylab=ylab, xlab="chromosomes", ylim=ylimit, xaxt="n", xaxs="i")
        abline(h=0)
        for (j in 2:max(chrom))
            abline(v=chrom.ends[j-1], lty='dashed')
        ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
        axis(side=1, at=ax, labels=chrom.labels, cex=.2, lwd=.5, las=1, cex.axis=1, cex.lab=1)
        amps <- data[,3]
        amps[amps>=5] <- 5.15
        amps[amps<5] <- NA
        points(pos, amps, pch=24, col='blue', bg='blue', cex=0.5)
        dels <- data[,3]
        dels[dels<=-2] <- -2.15
        dels[dels>-2] <- NA
        points(pos, dels, pch=25, col='red', bg='red', cex=0.5)
        ### MAD
        mad.value <- round(mad(copynumber(x)[chromosomes(x) < 23,i], na.rm=TRUE), digits=2)
        mtext(paste('MAD =', mad.value), side=3, line=0, adj=1)
        ### number of data points
        str <- paste(round(nclone / 1000), 'k x ', sep='')
        probe <- median(bpend(x)-bpstart(x)+1)
        if (probe < 1000) {
            str <- paste(str, probe, ' bp', sep='')
        } else {
            str <- paste(str, round(probe / 1000), ' kbp', sep='')
        }
        mtext(str, side=3, line=0, adj=0)
    }
})

setMethod("plot", signature(x="cghSeg", y="missing"),
function (x, y, dotres=1, ylimit=c(-2,5), ylab=expression(log[2]~ratio), build="GRCh37",... )
{
    for (i in 1:ncol(x)) {
        cat("Plotting sample", sampleNames(x)[i], "\n")
        segment         <- CGHbase:::.makeSegments(cbind(chromosomes(x), segmented(x)[,i]))
        chrom           <- chromosomes(x)
        data            <- data.frame(chrom, bpstart(x), copynumber(x)[,i])
        colnames(data)  <- c("chromosome", "position", "ratio")
        chrom.labels    <- as.character(unique(chrom))
        pos             <- bpstart(x)
        chrom.ends <- .getCumulativeChromosomeEnds(build)[1:max(chrom)]
        for (j in 2:max(chrom))
            pos[chrom == j] <- pos[chrom == j] + chrom.ends[j-1]
        nclone <- length(chrom)
        whichtoplot <- seq(1,nclone,by=dotres) #added 15/06/2009
        plot(pos[whichtoplot], data[whichtoplot,3], cex=.1, main=sampleNames(x)[i], ylab=ylab, xlab="chromosomes", ylim=ylimit, xaxt="n", xaxs="i")
        if (dotres != 1)
            mtext(paste('Plot resolution: ', 100/dotres, '%', sep=''), side=3, line=0)
        abline(h=0) 
        for (j in 2:max(chrom))
            abline(v=chrom.ends[j-1], lty='dashed')
        ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
        axis(side=1,at=ax,labels=chrom.labels,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
        for (jjj in (1:nrow(segment))) {
            segments(pos[segment[jjj,2]], segment[jjj,1], pos[segment[jjj,3]], segment[jjj,1], col="chocolate", lwd=3)        
        }
        amps <- data[,3]
        amps[amps>=5] <- 5.15
        amps[amps<5] <- NA
        points(pos, amps, pch=24, col='blue', bg='blue', cex=0.5)
        dels <- data[,3]
        dels[dels<=-2] <- -2.15
        dels[dels>-2] <- NA
        points(pos, dels, pch=25, col='red', bg='red', cex=0.5)
        ### MAD
        mad.value <- round(mad(copynumber(x)[chromosomes(x) < 23,i], na.rm=TRUE), digits=2)
        mtext(paste('MAD =', mad.value), side=3, line=0, adj=1)
        ### number of data points
        str <- paste(round(nclone / 1000), 'k x ', sep='')
        probe <- median(bpend(x)-bpstart(x)+1)
        if (probe < 1000) {
            str <- paste(str, probe, ' bp', sep='')
        } else {
            str <- paste(str, round(probe / 1000), ' kbp', sep='')
        }
        mtext(str, side=3, line=0, adj=0)
    }
})

setMethod("plot", signature(x="cghCall", y="missing"),
function (x, y, dotres=1, ylimit=c(-5,5), ylab=expression(log[2]~ratio), gaincol='blue', losscol='red', misscol=NA, build='GRCh37',... )
{
    calls           <- calls(x)
    nsamples        <- ncol(x)
    if (!is.null(probamp(x))) nclass <- 4 else nclass     <- 3
    chrom           <- chromosomes(x)
    pos             <- bpstart(x)
    pos2            <- bpend(x)
    chrom.ends <- .getCumulativeChromosomeEnds(build)[1:max(chrom)]
    for (j in 2:max(chrom)) {
        pos[chrom == j] <- pos[chrom == j] + chrom.ends[j-1]
        pos2[chrom == j] <- pos2[chrom == j] + chrom.ends[j-1]
    }
    nclone          <- length(chrom)

    for (i in 1:ncol(x)) {
    
        cat("Plotting sample", sampleNames(x)[i], "\n")
    
        genomdat        <- copynumber(x)[,i]
        probsdraw       <- cbind(probloss(x)[,i], probnorm(x)[,i], probgain(x)[,i])
        if (!is.null(probamp(x))) probsdraw <- cbind(probsdraw, probamp(x)[,i])
        lt              <- 0
        
        if (nclass==4) {
            ticksamp    <- which(probsdraw[,nclass] >= 0.5)
            lt          <- length(ticksamp)        
            probsdraw   <- cbind(probsdraw[,1:2], probsdraw[,3] + probsdraw[,4])
        }
        
        segment         <- CGHbase:::.makeSegments(cbind(chromosomes(x), segmented(x)[,i]))
        
        widths          <- segment[,3] - segment[,2] + 1
        par(mar=c(5, 4, 4, 4) + 0.2)
        
        ### Plot the probability bars
        plot(NA, xlim=c(0, max(pos2)), ylim=c(0,1), xlab=NA, ylab=NA, las=1, xaxs='i', xaxt='n', yaxs='i')
        if (!is.na(misscol)) {
            rect(0, -1, max(pos2), 1, col=misscol, border=NA)
            rect(pos, -1, pos2, 1, col='white', border=NA)
        }
        rect(pos[segment[,2]], 0, pos2[segment[,3]], probsdraw[segment[,2],1], col=losscol, border=losscol)
        rect(pos[segment[,2]], 1, pos2[segment[,3]], 1-probsdraw[segment[,2],3], col=gaincol, border=gaincol)
        
        lim <- par("usr")
        lim[3:4] <- ylimit
        par(usr=lim)
        dticks <- seq(ylimit[1], ylimit[2], by=1)
        axis(4, at=dticks, labels=dticks, srt=270, las=1, cex.axis=1, cex.lab=1)
        if (lt > 0) {
            axis(3,at=pos[ticksamp], labels=FALSE, col=gaincol, col.axis="black", srt=270, las=1, cex.axis=1, cex.lab=1)
        }     
        box()   
        abline(h=0)
        
        ### Add axis labels
        mtext(ylab, side=4, line=3, srt=270)
        mtext("probability", side=2, line=3, srt=270)
        mtext('chromosomes', side=1, line=3)
        
        #### add vert lines at chromosome ends
        for (j in 2:max(chrom))
            abline(v=chrom.ends[j-1], lty='dashed')

        title(sampleNames(x)[i])
        if (dotres != 1)
            mtext(paste('Plot resolution: ', 100/dotres, '%', sep=''), side=3, line=0)
        
        ### Add log2ratios
        if (dotres>0) {
            whichtoplot <- seq(1,nclone,by=dotres) #added 15/06/2009
            points(pos[whichtoplot],genomdat[whichtoplot],cex=.1)
        }
        
        ### X-axis with chromosome labels
        ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
        axis(side=1,at=ax,labels=unique(chrom),cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
            
        ### segment means
        for (jjj in (1:nrow(segment)))
            segments(pos[segment[jjj,2]], segment[jjj,1], pos[segment[jjj,3]], segment[jjj,1], col="chocolate", lwd=3)        

        ### MAD
        mad.value <- round(mad(copynumber(x)[chromosomes(x) < 23,i], na.rm=TRUE), digits=2)
        mtext(paste('MAD =', mad.value), side=3, line=0, adj=1)

        ### number of data points
        str <- paste(round(nclone / 1000), 'k x ', sep='')
        probe <- median(bpend(cgh)-bpstart(cgh)+1)
        if (probe < 1000) {
            str <- paste(str, probe, ' bp', sep='')
        } else {
            str <- paste(str, round(probe / 1000), ' kbp', sep='')
        }
        mtext(str, side=3, line=0, adj=0)
    }
})

setMethod("frequencyPlot", signature(x="cghCall", y="missing"),
function (x, y, main='Frequency Plot', gaincol='blue', losscol='red', misscol=NA, build='GRCh37',... )
{
  chrom <- chromosomes(x)
  pos <- bpstart(x)
  pos2 <- bpend(x)
  chrom.ends <- .getCumulativeChromosomeEnds(build)[1:max(chrom)]
  for (j in 2:max(chrom)) {
    pos[chrom == j] <- pos[chrom == j] + chrom.ends[j-1]
    pos2[chrom == j] <- pos2[chrom == j] + chrom.ends[j-1]
  }
  calls <- calls(x)
  loss.freq <- rowMeans(calls < 0)
  gain.freq <- rowMeans(calls > 0)
  plot(NA, xlim=c(0, max(pos2)), ylim=c(-1,1), type='n', xlab='chromosomes', ylab='frequency', xaxs='i', xaxt='n', yaxs='i', yaxt='n', main=main,...)
  if (!is.na(misscol)) {
    rect(0, -1, max(pos2), 1, col=misscol, border=NA)
    rect(pos, -1, pos2, 1, col='white', border=NA)
  }
  rect(pos, 0, pos2, gain.freq, col=gaincol, border=gaincol)
  rect(pos, 0, pos2, -loss.freq, col=losscol, border=losscol)
  box()
  abline(h=0)
  for (j in 2:max(chrom))
    abline(v=chrom.ends[j-1], lty='dashed')
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1,at=ax,labels=unique(chrom),cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
  mtext('gains', side=2, line=3, at=0.5)
  mtext('losses', side=2, line=3, at=-0.5)
  ### number of data points
  str <- paste(round(nrow(x) / 1000), 'k x ', sep='')
  probe <- median(bpend(x)-bpstart(x)+1)
  if (probe < 1000) {
    str <- paste(str, probe, ' bp', sep='')
  } else {
    str <- paste(str, round(probe / 1000), ' kbp', sep='')
  }
  mtext(str, side=3, line=0, adj=0)
})

setMethod("frequencyPlot", signature(x="cghRegions", y="missing"),
function (x, y, main='Frequency Plot', gaincol='blue', losscol='red', misscol=NA, build='GRCh37',... )
{
  chrom <- chromosomes(x)
  pos <- bpstart(x)
  pos2 <- bpend(x)
  chrom.ends <- .getCumulativeChromosomeEnds(build)[1:max(chrom)]
  for (j in 2:max(chrom)) {
    pos[chrom == j] <- pos[chrom == j] + chrom.ends[j-1]
    pos2[chrom == j] <- pos2[chrom == j] + chrom.ends[j-1]
  }
  calls <- regions(x)
  loss.freq <- rowMeans(calls < 0)
  gain.freq <- rowMeans(calls > 0)
  plot(NA, xlim=c(0, max(pos2)), ylim=c(-1,1), type='n', xlab='chromosomes', ylab='frequency', xaxs='i', xaxt='n', yaxs='i', yaxt='n', main=main,...)
  if (!is.na(misscol)) {
    rect(0, -1, max(pos2), 1, col=misscol, border=NA)
    rect(pos, -1, pos2, 1, col='white', border=NA)
  }
  rect(pos, 0, pos2, gain.freq, col=gaincol, border=gaincol)
  rect(pos, 0, pos2, -loss.freq, col=losscol, border=losscol)
  box()
  abline(h=0)
  for (j in 2:max(chrom))
    abline(v=chrom.ends[j-1], lty='dashed')
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1,at=ax,labels=unique(chrom),cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
  mtext('gains', side=2, line=3, at=0.5)
  mtext('losses', side=2, line=3, at=-0.5)
  ### number of data points
  mtext(paste(nrow(x), 'regions'), side=3, line=0, adj=0)
})

make_cghRawPlus <- function(input) {
    if (class(input) == "character") input  <- read.table(input, header=T, sep="\t", fill=T, quote="")
    if (class(input[,2]) == 'factor')
        input[,2] <- as.character(input[,2])
    if (class(input[,2]) == 'character') {
        input[,2] <- sub('^chr', '', input[,2])
        input[input[,2] == 'X', 2] <- '23'
        input[input[,2] == 'Y', 2] <- '24'
        input[input[,2] == 'MT', 2] <- '25'
        input[,2] <- as.integer(input[,2])
    }
    input <- input[order(input[,2], input[,3]),]
    copynumber  <- as.matrix(input[,5:ncol(input)])
    rownames(copynumber) <- input[,1]
    if (ncol(copynumber) == 1)
      colnames(copynumber) <- colnames(input)[5]
    annotation  <- data.frame(Chromosome=input[,2], Start=input[,3], End=input[,4], row.names=input[,1])
    metadata    <- data.frame(labelDescription=c("Chromosomal position", "Basepair position start", "Basepair position end"), row.names=c("Chromosome", "Start", "End"))    
    dimLabels   <- c("featureNames", "featureColumns")
    annotation  <- new("AnnotatedDataFrame", data=annotation, dimLabels=dimLabels, varMetadata=metadata)   
    result      <- new("cghRaw", copynumber=copynumber, featureData=annotation)
}
environment(make_cghRawPlus) <- environment(CGHcall:::make_cghRaw)
make_cghRaw <- make_cghRawPlus

preprocessPlus <- 
function (input, maxmiss = 30, nchrom = 23, ...) 
{
    input <- input[!is.na(chromosomes(input)) & !is.na(bpstart(input))]
    input <- input[chromosomes(input) != 0 & bpstart(input) != 
        0, ]
    input <- input[chromosomes(input) <= nchrom, ]
    countMissing <- function(x, maxMissing) {
        if (length(x[is.na(x)]) <= maxMissing) 
            return(TRUE)
        return(FALSE)
    }
    nummis <- length(copynumber(input)[is.na(copynumber(input))])
    if (nummis > 0) {
        allowed <- ncol(input) * maxmiss/100
        filter <- apply(copynumber(input), 1, countMissing, allowed)
        input <- input[filter, ]
    }
    nummis <- length(copynumber(input)[is.na(copynumber(input))])
    if (nummis > 0) {
        if (!exists("k")) 
            k <- 10
        new.k <- ceiling((1 - maxmiss/100) * ncol(input))
        new.k <- min(new.k, k)
        if (new.k != 10) 
            cat("Changing impute.knn parameter k from", k, "to", 
                new.k, "due to small sample size.\n")
        if (new.k > 1) {
            nrows <- nrow(copynumber(input))
            nrpi <- floor(nrows/100)
            resultsall <- c()
            for (i in 1:(100 - 1)) {
                datforimp <- copynumber(input)[((i - 1) * nrpi + 
                  1):(nrpi * i), ]
                result <- impute::impute.knn(datforimp, k = new.k)
                resultsall <- rbind(resultsall, result$data)
            }
            datforimp <- copynumber(input)[((100 - 1) * nrpi + 
                1):nrows, ]
            result <- impute::impute.knn(datforimp, k = new.k)
            resultsall <- rbind(resultsall, result$data)
            copynumber(input) <- CGHcall:::.assignNames(resultsall, 
                input)
        }
    }
    input
}
environment(preprocessPlus) <- environment(CGHcall:::preprocess)
preprocess <- preprocessPlus

.getCumulativeChromosomeEnds <- function(build) {
    build <- as.integer(gsub('[^0-9]', '', build))
    if (build == 34 || build == 16) {
       chromosome.sizes <- c(246127941, 243615958, 199344050, 191731959, 181034922, 170914576, 158545518, 146308819, 136372045, 135037215, 134482954, 132078379, 113042980, 105311216, 100256656, 90041932, 81860266, 76115139, 63811651, 63741868, 46976097, 49396972, 153692391, 50286555)
    } else if (build == 35 || build == 17) {
       chromosome.sizes <- c(245522847, 243018229, 199505740, 191411218, 180857866, 170975699, 158628139, 146274826, 138429268, 135413628, 134452384, 132449811, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49554710, 154824264, 57701691)
    } else if (build == 36 || build == 18) {
       chromosome.sizes <- c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
    } else {
       chromosome.sizes <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
    }
    for (i in 2:length(chromosome.sizes))
      chromosome.sizes[i] <- chromosome.sizes[i] + chromosome.sizes[i-1]
    chromosome.sizes
}

.getCentromerePlus <- function(build) {
    build <- as.integer(gsub('[^0-9]', '', build))
    centromere <- matrix(NA, 23, 2);
    if (build == 34 || build == 16) {
        centromere[1,] <- c(120093537, 122333537)
        centromere[2,] <- c(91810927, 94810927)
        centromere[3,] <- c(90425755, 93325755)
        centromere[4,] <- c(49575659, 52575659)
        centromere[5,] <- c(46451142, 49451142)
        centromere[6,] <- c(58877002, 61877002)
        centromere[7,] <- c(57832528, 60832528)
        centromere[8,] <- c(43856263, 46856263)
        centromere[9,] <- c(44443361, 47443361)
        centromere[10,] <- c(39208941, 41588941)
        centromere[11,] <- c(51602318, 54602318)
        centromere[12,] <- c(34747961, 36142961)
        centromere[13,] <- c(14900000, 16768000)
        centromere[14,] <- c(15070000, 18070000)
        centromere[15,] <- c(15260000, 18260000)
        centromere[16,] <- c(36366889, 38166889)
        centromere[17,] <- c(22408570, 25408570)
        centromere[18,] <- c(15398887, 16762885)
        centromere[19,] <- c(26923622, 29923622)
        centromere[20,] <- c(26314569, 29314569)
        centromere[21,] <- c(10260000, 13260000)
        centromere[22,] <- c(11330000, 14330000)
        centromere[23,] <- c(57548803, 60548803)
        # centromere[24,] <- c(9757849, 12757849)
    } else if (build == 35 || build == 17) {
        centromere[1,] <- c(121147476, 123387476)
        centromere[2,] <- c(91748045, 94748045)
        centromere[3,] <- c(90587544, 93487544)
        centromere[4,] <- c(49501045, 52501045)
        centromere[5,] <- c(46441398, 49441398)
        centromere[6,] <- c(58938125, 61938125)
        centromere[7,] <- c(57864988, 60864988)
        centromere[8,] <- c(43958052, 46958052)
        centromere[9,] <- c(46035928, 49035928)
        centromere[10,] <- c(39244941, 41624941)
        centromere[11,] <- c(51450781, 54450781)
        centromere[12,] <- c(34747961, 36142961)
        centromere[13,] <- c(16000000, 17868000)
        centromere[14,] <- c(15070000, 18070000)
        centromere[15,] <- c(15260000, 18260000)
        centromere[16,] <- c(35143302, 36943302)
        centromere[17,] <- c(22187133, 22287133)
        centromere[18,] <- c(15400898, 16764896)
        centromere[19,] <- c(26923622, 29923622)
        centromere[20,] <- c(26267569, 28033230)
        centromere[21,] <- c(10260000, 13260000)
        centromere[22,] <- c(11330000, 14330000)
        centromere[23,] <- c(58465033, 61465033)
        # centromere[24,] <- c(11237315, 12237315)
    } else if (build == 36 || build == 18) {
        centromere[1,] <- c(121236957, 123476957)
        centromere[2,] <- c(91689898, 94689898)
        centromere[3,] <- c(90587544, 93487544)
        centromere[4,] <- c(49354874, 52354874)
        centromere[5,] <- c(46441398, 49441398)
        centromere[6,] <- c(58938125, 61938125)
        centromere[7,] <- c(58058273, 61058273)
        centromere[8,] <- c(43958052, 46958052)
        centromere[9,] <- c(47107499, 50107499)
        centromere[10,] <- c(39244941, 41624941)
        centromere[11,] <- c(51450781, 54450781)
        centromere[12,] <- c(34747961, 36142961)
        centromere[13,] <- c(16000000, 17868000)
        centromere[14,] <- c(15070000, 18070000)
        centromere[15,] <- c(15260000, 18260000)
        centromere[16,] <- c(35143302, 36943302)
        centromere[17,] <- c(22187133, 22287133)
        centromere[18,] <- c(15400898, 16764896)
        centromere[19,] <- c(26923622, 29923622)
        centromere[20,] <- c(26267569, 28033230)
        centromere[21,] <- c(10260000, 13260000)
        centromere[22,] <- c(11330000, 14330000)
        centromere[23,] <- c(58598737, 61598737)
    } else { # 37 / 19
        centromere[1,] <- c(121535434, 124535434)
        centromere[2,] <- c(92326171, 95326171)
        centromere[3,] <- c(90504854, 93504854)
        centromere[4,] <- c(49660117, 52660117)
        centromere[5,] <- c(46405641, 49405641)
        centromere[6,] <- c(58830166, 61830166)
        centromere[7,] <- c(58054331, 61054331)
        centromere[8,] <- c(43838887, 46838887)
        centromere[9,] <- c(47367679, 50367679)
        centromere[10,] <- c(39254935, 42254935)
        centromere[11,] <- c(51644205, 54644205)
        centromere[12,] <- c(34856694, 37856694)
        centromere[13,] <- c(16000000, 19000000)
        centromere[14,] <- c(16000000, 19000000)
        centromere[15,] <- c(17000000, 20000000)
        centromere[16,] <- c(35335801, 38335801)
        centromere[17,] <- c(22263006, 25263006)
        centromere[18,] <- c(15460898, 18460898)
        centromere[19,] <- c(24681782, 27681782)
        centromere[20,] <- c(26369569, 29369569)
        centromere[21,] <- c(11288129, 14288129)
        centromere[22,] <- c(13000000, 16000000)
        centromere[23,] <- c(58632012, 61632012)
        # centromere[24,] <- c(10104553, 13104553)
    }
    centromere <- apply(centromere, 1, mean);
    return(centromere);
}

.convertChromosomeToArmPlus <- function(dataframe, build) { #changed 22/06/2009; 
    cat("Dividing chromosomes into arms using centromere positions from", build, "\n\n");
    centromere  <- .getCentromerePlus(build);
    chr <- dataframe[,2]
    bp <- dataframe[,3]
    chrlev <- unique(chr)
    a<-1
    chrarms <- c()
    for(i in chrlev){
    print(i)
    chri <- which(chr==i)
    bpi <- bp[chri]
    wbpi <- length(which(bpi<=centromere[i]))
    wbpil <- length(which(bpi>centromere[i]))
    if(wbpi>0) {chrarms <- c(chrarms,rep(a,wbpi));a<-a+1}
    if(wbpil>0) {chrarms <- c(chrarms,rep(a,wbpil));a<-a+1}  
    }  
    dataframe[,2] <- chrarms
    return(dataframe)
}

CGHcallPlus <- function(inputSegmented, prior="auto", nclass=4, organism="human", robustsig="yes",nsegfit=3000,maxnumseg=100,minlsforfit=0.5, build="GRCh37") {
  #  library(CGHcall)
#    setwd("C:\\VUData\\Maria")
#    #rm(list=ls())
#    memory.limit(size = 4000)
#    load("seg.Rdata")
#    inputSegmented <- seg; prior="auto"; nclass=4; organism="human"; robustsig="yes";nsegfit=500;maxnumseg=100;minlsforfit=0.5
    timeStarted <- proc.time()
    print("changed")
    gc() #memory usage in Mb
    
 
    ## Author: Mark van de Wiel
    ## Maintainer: Mark van de Wiel
    ## Email: mark.vdwiel@vumc.nl
 
    
    ## Changes since previous version
    ## - input is of class cgh
   
    timeStarted <- proc.time()

    ## Read input data
    #normalizedData  <- copynumber(inputSegmented)
    #combined        <- data.frame(featureNames(inputSegmented), chromosomes(inputSegmented), bpstart(inputSegmented), copynumber(inputSegmented),
    #segmented(inputSegmented))
    bw <- round(20*minlsforfit*nrow(copynumber(inputSegmented))/44000) 
    ncolscl     <- ncol(copynumber(inputSegmented))
  #  datareg     <- CGHcall:::.MakeData(data.frame(featureNames(inputSegmented), chromosomes(inputSegmented), bpstart(inputSegmented), copynumber(inputSegmented), segmented(inputSegmented)))  
    whna <- !is.na(chromosomes(inputSegmented)) & !is.na(featureNames(inputSegmented))
    datmat <- cbind(copynumber(inputSegmented), segmented(inputSegmented))[whna,]
    #datmat <- cbind(c(-1.2,2.1,-0.3,-0.4),c(0,0,0,0));chr<-c(1,1,1,1);nc<-1
    posit       <- bpstart(inputSegmented)[whna]
    chr         <- chromosomes(inputSegmented)[whna]
    naam        <- featureNames(inputSegmented)[whna]
    rm(whna,inputSegmented);gc()
    nc    <- ncol(datmat)/2 
    nclone <-nrow(datmat)
    datareg  <- .MakeData(datmat[,-(1:nc),drop=FALSE],chr)  
    #datareg  <- CGHcall:::.MakeData(datmat[,-(1:nc),drop=FALSE],chr)  
   
    #rm(combined); 
    
    gc() #added 24/11/2009
              
    #datareg <- .MakeData(combined)
    #save(datareg,file="datareg.Rdata")
    #load("datareg.Rdata")
     dataprob    <- data.frame(naam, chr, posit)
    
 
    
    ## Determine method for prior probabilities
    if (prior == 'auto') {
        if (nc >  20) prior <- "not all";
        if (nc <= 20) prior <- "all";
    }
    
    ### Convert to chromosome arm if neccessary
    if (prior != "all" && organism == "human") {
       # temp            <- datareg[[1]];
        dataprob    <- .convertChromosomeToArmPlus(dataprob, build); 
       # datareg[[1]]    <- .convertChromosomeToArmPlus(temp, build);
        chr <- dataprob[,2] #BUG;repaired 22/06/09
    }

    cat("EM algorithm started ... \n");
    
    dat         <- c()
    datsm       <- c()
    datall         <- c()
    datsmall       <- c()
    regions     <- c()
    regionsall  <- c()
    regionsdat  <- c()
    regionsdatall  <- c()
    nclones     <- length(chr)
     profile     <- c()
     profileall <- c()
    countreg    <- 0
    newregions <- list()
#    nclonesj <- 0
    thr <- 2
   
    #here is the selection of regions used for fitting. First select regions long enough and with smoothed signal <= thr. Then, maximize 
    #regions per prof to maxseg
    #nc<-2
    gc() 
    datall <- as.vector(datmat[,1:nc]) #changed 16/7/10
    datsmall      <- as.vector(datmat[,-(1:nc)])
    rm(datmat);gc();
  
    for (j in (1:nc)) {
        regions1    <- datareg[[j]]
        nreg1       <- nrow(regions1)
        profileall     <- c(profileall, rep(j,nreg1))
        regionsall     <- rbind(regionsall, regions1)
        toaddall       <- (j-1)*c(nclones, nclones)
        regionsdatall  <- rbind(regionsdatall, regions1[,-3]+toaddall)

        segl        <- regions1[,2]-regions1[,1]+1
        whseg       <- which(segl>=bw)
        regions2    <- regions1[whseg,,drop=FALSE] #selects regions long enough to participate in fitting
        regions2    <- regions2[which(abs(regions2[,3]) <= thr),,drop=FALSE] #selects regions with smoothed signal small enough
        nreg2       <- nrow(regions2)
        countreg <- countreg + min(maxnumseg,nreg2)
        newregions[[j]] <- regions2
        }
        fracforfit <- min(1,nsegfit/countreg)
        
        
        for (j in (1:nc)) {
        regions2 <- newregions[[j]]
        nreg2       <- nrow(regions2)
        ordsig      <- order(regions2[,3])
        regions2b   <- data.frame(regions2,genord = 1:nreg2)
        regions2c   <- regions2b[ordsig,,drop=FALSE]
        if(nreg2 >= maxnumseg | fracforfit<1){
            takeseg     <- floor(seq(1,nreg2,length.out=round(min(nreg2,maxnumseg)*fracforfit)))
            regions2c <- regions2c[takeseg,,drop=FALSE]
            }
        regions2d <- regions2c[order(regions2c[,4]),,drop=FALSE]
        
        nreg2       <- nrow(regions2d)
        profile    <- c(profile, rep(j,nreg2))
        regions     <- rbind(regions, regions2d)
        #takerow <- c();for(k in 1:nrow(regions2d)) {takerow <- c(takerow,regions2d[k,1]:regions2d[k,2])}
#        dat <- c(dat,datareg[[1]][takerow,(3+j)])
#        datsm <- c(datsm,datareg[[1]][takerow,(3+nc+j)])
        toadd       <- (j-1)*c(nclones, nclones)
        regionsdat  <- rbind(regionsdat, regions2d[,-(3:4)]+toadd) #delete smrat value and index
        #        nclonesj <- sum((regions2d[,2]+1)-regions2d[,1]) 
    }
    gc()

    takechr     <- function(reg, chrom) {
        chrom[reg[1]]
    }
    
    chrarmreg   <- as.vector(apply(as.matrix(regions), 1, takechr, chrom=chr))
    
    nreg        <- nrow(regionsdat)
    nregall     <- nrow(regionsdatall)
    print(paste("Total number of segments present in the data:",nregall))
    print(paste("Number of segments used for fitting the model:",nreg))
    allnc       <- sapply(1:nreg, CGHcall:::.countcl, regionsdat=regionsdat)
    allsum      <- sapply(1:nreg, CGHcall:::.sumreg, dat=datall, regionsdat=regionsdat)
    allsumsq    <- sapply(1:nreg, CGHcall:::.sumsqreg, dat=datall, regionsdat=regionsdat)
    varregall   <- sapply(1:nreg, CGHcall:::.varregtimescount, counts=allnc, dat=datall, regionsdat=regionsdat)
    lev         <- 1:nc
    varprofnc   <- cbind(varregall, profile, allnc)
    varprof     <- sapply(lev, CGHcall:::.varproffun, vcnmat=varprofnc, profile=profile)
    
    gc()

    selk <- function(k, varprof,prof) {  #changed 22/06/09
        profk <- prof[k]
        return(max(0.001,varprof[profk]))
    }
    varprofall <- sapply(1:nreg, selk, varprof=varprof,prof=profile)

    allmean     <- allsum/allnc
    allmeanl    <- allmean[allmean < -0.15 & allmean > -0.7]
    allmeang    <- allmean[allmean>0.15 & allmean<0.7]
    allmean0    <- allmean[abs(allmean) <= 0.15]
    sdlst       <- if(length(allmeanl) <= 1) 0.01 else mad(allmeanl)
    meanlst     <- if(length(allmeanl) <= 0) -0.3 else mean(allmeanl)
    sdgst       <- if(length(allmeang) <= 1) 0.01 else mad(allmeang)
    meangst     <- if(length(allmeang) <= 0) 0.3 else mean(allmeang)
    sd0st       <- if(length(allmean0) <= 1) 0.01 else mad(allmean0)
    #sd0st       <- if(length(allmean0) <= 1) 0.01 else sd(allmean0)

    profchrom   <- chrarmreg
    if(robustsig=="no") {
        bstart      <- c(-log(0.2), -log(-(meanlst+0.1)), sqrt(-log(0.5)), -log(meangst-0.1), -log(0.2), 0.0001, sdlst, sd0st, sdgst, 0.0001, 0.0001)
    } else {
        bstart      <- c(-log(0.2), -log(-(meanlst+0.1)), sqrt(-log(0.5)), -log(meangst-0.1), -log(0.2), 0.0001, sdlst, sd0st+max(0,sqrt(sdgst^2/8+sdlst^2/8)-sd0st), sdgst, 0.0001, 0.0001) #robust option 
    }
    mus         <- c(-0.10-exp(-bstart[2])-0.3 - exp(-bstart[1]), -0.10-exp(-bstart[2]), -0.05+0.1*exp(-(bstart[3])^2), 0.1 + exp(-bstart[4]), 2*(0.10 + exp(-bstart[4]))+0.3+exp(-bstart[5]))
    priorp      <- matrix(rep(c(0.01, 0.09, 0.8, 0.08, 0.01, 0.01), nreg), ncol=6, byrow=TRUE)
    alpha0      <- CGHcall:::.alpha0all(nreg, profchrom, priorp, bstart, varprofall, allsum, allsumsq, allnc,robustsig,prior)

    maxiter     <- 10
    stop        <- 0
    iter        <- 1
    thrpara     <- 0.015 #changed 15/6/09
    
    gc()
    
    while (stop == 0 & iter <= maxiter) {
        print(gc())
        cat("Calling iteration", iter, ":\n")
        posterior0  <- sapply(1:nreg, CGHcall:::.posteriorp, priorp=alpha0, pm=bstart, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig)
        likprev     <- CGHcall:::.totallik(bstart, nreg=nreg, posteriorprev=posterior0, alphaprev=alpha0, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig)
        nlmres      <- nlm(CGHcall:::.totallik, bstart, nreg=nreg, posteriorprev=posterior0, alphaprev=alpha0, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig)
        bprev       <- bstart
        bstart      <- nlmres$est
        alpha0      <- CGHcall:::.alpha0all(nreg, profchrom, alpha0, bstart, varprofall, allsum, allsumsq, allnc, robustsig, prior)
        rl          <- CGHcall:::.reallik4(nreg, alpha0, bstart, varprofall, allsum, allsumsq, allnc, robustsig)
        llmin       <- nlmres$min
        musprev     <- c(-0.10-exp(-bprev[2])-0.3 - exp(-bprev[1]), -0.10-exp(-bprev[2]), -0.05+0.1*exp(-(bprev[3])^2), 0.1 + exp(-bprev[4]),(log2(2)/log2(1.5))*(0.1 + exp(-bprev[4])), (log2(2)/log2(1.5))*(0.10 + exp(-bprev[4]))+0.3+exp(-bprev[5]))
        mus         <- c(-0.10-exp(-bstart[2])-0.3 - exp(-bstart[1]), -0.10-exp(-bstart[2]), -0.05+0.1*exp(-(bstart[3])^2), 0.1 + exp(-bstart[4]),(log2(2)/log2(1.5))*(0.1 + exp(-bstart[4])),(log2(2)/log2(1.5))*(0.10 + exp(-bstart[4]))+0.3+exp(-bstart[5]))
        param       <- c(mus, bstart[6], bstart[7], bstart[8], bstart[9], bstart[10],bstart[11])
        paramprev   <- c(musprev, bprev[6], bprev[7], bprev[8], bprev[9], bprev[10],bprev[11])
     
        if(robustsig=="no") {
            printmat <- matrix(c(j,rl,mus,sqrt(bstart[7]^2+(bstart[6])^2 + 0.0001),bstart[7],bstart[8],bstart[9],sqrt((bstart[10])^2+(bstart[9])^2 + 0.0001),sqrt((bstart[11])^2 + (bstart[10])^2+ (bstart[9])^2+ 0.0001)),nrow=1)
            colnames(printmat) <- c("j","rl","mudl","musl","mun","mug","mudg","mua","sddl","sdsl","sdn","sdg","sddg","sda")
            print(printmat)
        } else {
            printmat<-matrix(c(j,rl,mus,sqrt(bstart[7]^2+(bstart[6])^2 + 0.0001),bstart[7],sqrt(1/8*((bstart[7])^2 + 0.0001)+1/8*((bstart[9])^2 + 0.0001)+((bstart[8]^2)+0.0001)),bstart[9],sqrt((bstart[10])^2+(bstart[9])^2 + 0.0001),sqrt((bstart[11])^2 + (bstart[10])^2 + (bstart[9])^2 + 0.0001)),nrow=1)  #robust option
            colnames(printmat) <- c("j","rl","mudl","musl","mun","mug","mudg","mua","sddl","sdsl","sdn","sdg","sddg","sda")
            print(printmat)
        }
        if (max(abs(param[c(2:4,8)]-paramprev[c(2:4,8)])) <= thrpara) stop <- 1 #changed 15/6/09; also sdn should converge
        iter <- iter+1
    }
    mudl<-printmat[3];mul <- printmat[4];mug<-printmat[6];mudg<- printmat[7];mua<-printmat[8];sdl<-printmat[10];sdg<-printmat[12] #added 17/12
   # pm <- as.vector(printmat)
#    plot(function(x) dnorm(x,mean=pm[4],sd=pm[11]),-1,1,col="red")
#   plot(function(x) dnorm(x,mean=pm[5],sd=pm[12]),-1,1,col="black",add=TRUE)
#   plot(function(x) dnorm(x,mean=pm[6],sd=pm[13]),-1,1,col="green",add=TRUE)
    cat("EM algorithm done ...\n")
    best            <- bstart
    
    #now start computing posteriors for ALL data
    cat("Computing posterior probabilities for all segments ...\n")
    #nregall        <- nrow(regionsdatall)
    allncall       <- sapply(1:nregall, CGHcall:::.countcl, regionsdat=regionsdatall)
    allsumall      <- sapply(1:nregall, CGHcall:::.sumreg, dat=datall, regionsdat=regionsdatall)
    allmeanall     <- allsumall/allncall  
    allsumsqall    <- sapply(1:nregall, CGHcall:::.sumsqreg, dat=datall, regionsdat=regionsdatall)
    varregallall   <- sapply(1:nregall, CGHcall:::.varregtimescount, counts=allncall, dat=datall, regionsdat=regionsdatall)
    lev             <- 1:nc
    varprofncall   <- cbind(varregallall, profileall, allncall)
    varprof_all     <- sapply(lev, CGHcall:::.varproffun, vcnmat=varprofncall, profile=profileall)
    varprof_allall <- sapply(1:nregall, selk, varprof=varprof_all,prof=profileall)
    gc()
    
    if(prior!="all"){
    chrarm_alpha0 <- cbind(chrarmreg,alpha0)
    chrarm_alpha0uni <- unique(chrarm_alpha0) 
    chrarmregall   <- as.vector(apply(as.matrix(regionsall), 1, takechr, chrom=chr))
    alpha0_all      <- t(sapply(chrarmregall,function(x){chrarm_alpha0uni[chrarm_alpha0uni[,1]==x,-1]})) 
    } else alpha0_all <- matrix(rep(alpha0[1,],nregall),byrow=T,nrow=nregall)
    
    posteriorfin    <- t(sapply(1:nregall, CGHcall:::.posteriorp, priorp=alpha0_all, pm=best, varprofall=varprof_allall, allsum=allsumall, allsumsq=allsumsqall, allnc=allncall, robustsig=robustsig))
    #posteriorfin1   <- cbind(allncall, allmeanall, posteriorfin)  #REMOVED 17/12
    #amps            <- posteriorfin1[posteriorfin1[,5]+posteriorfin1[,6]>=0.5,] #REMOVED 17/12
    
    overrule <- function(id,mnall){
    #id <- 393;mnall<-allmeanall
    toret <- posteriorfin[id,]   
    mn <- mnall[id]
    row1 <- toret
    if((mn > mug+sdg) & (which.max(row1)<=3)) {
        sumr <- sum(row1[4],row1[5],row1[6]); #added 17/12/10; if max prob belongs to non-gain state, overrule; redistribute probability
        if(sumr > 10^(-10)) toret<- c(0,0,0,row1[4]/sumr,row1[5]/sumr,row1[6]/sumr) else {
            toret<- c(0,0,0,0,0,0)
            wm<-which.min(c(abs(mn-mug),abs(mn-mudg),abs(mn-mua)))+3
            toret[wm]<-1 #assigns probability 1 to status with mode closest to the mean segment value; only used if probability calculation fails.
        }
    }
    if((mn < mul-sdl) & (which.max(row1)>=3)) {
        sumr <- sum(row1[1],row1[2]); #added 17/12/10; if max prob belongs to non-loss state, overrule
         if(sumr > 10^(-10)) toret<- c(row1[1]/sumr,row1[2]/sumr,0,0,0,0) else {
            toret<- c(0,0,0,0,0,0)
            wm<-which.min(c(abs(mn-mudl),abs(mn-mul)))
            toret[wm]<-1 #assigns probability 1 to status with mode closest to the mean segment value; only used if probability calculation fails.
        } 
    }
    return(toret)
    }
    
    posteriorfin <- t(sapply(1:nregall,overrule,mnall=allmeanall)) #OVERRULE FUNCTION ADD 17/12/10 TO DEAL WITH EXTREMES
    
    amporgain <- function(row) {
        if (row[6] >= 0.5) return(c(row[1]+row[2], row[3], row[4]+row[5], row[6]))
        else return(c(row[1]+row[2], row[3], row[4]+row[5], row[6]))
    }
   
    if (nclass == 3) posteriorfin1 <- cbind(posteriorfin[,1]+posteriorfin[,2], posteriorfin[,3], posteriorfin[,4]+posteriorfin[,5]+posteriorfin[,6], posteriorfin[,6])
    if (nclass == 4) posteriorfin1 <- t(apply(posteriorfin, 1, amporgain))
    
   
    
    

    posteriorfin2   <- cbind(profileall, posteriorfin1)
    rm(posteriorfin,posteriorfin1);gc() #added24/11/2009
    regionsprof     <- cbind(profileall, regionsall)
    listcall <- list(posteriorfin2,nclone,nc,nclass,regionsprof,params = printmat[,-c(1,2)])
    #save(listcall,file="listcall.Rdata")
    gc()
    timeFinished <- round((proc.time() - timeStarted)[1] / 60)
    cat("Total time:", timeFinished, "minutes\n")
    return(listcall)
}
environment(CGHcallPlus) <- environment(CGHcall:::CGHcall)
CGHcall <- CGHcallPlus

CGHregionsPlus <- function(input, ...) {
  regions <- CGHregions:::CGHregions(input, ...)
  # End positions of regions should be the end position of the last data point of that region,
  # but instead CGHregions returns the start position of the last data point.
  # Check if that is indeed the case:
  if (class(input) == 'cghCall') {
    if (sum(regions@featureData@data$End %in% input@featureData@data$Start) > sum(regions@featureData@data$End %in% input@featureData@data$End))
      for (row in rownames(regions@featureData@data))
        regions@featureData@data[row, 'End'] <- input@featureData@data[input@featureData@data$Chromosome == regions@featureData@data[row, 'Chromosome'] & input@featureData@data$Start == regions@featureData@data[row, 'End'], 'End'][1]
  }
  regions
}
environment(CGHregionsPlus) <- environment(CGHregions:::CGHregions)
CGHregions <- CGHregionsPlus

# EOF
