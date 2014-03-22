# TOOL acgh-call-aberrations.R: "Call aberrations from segmented copy number data" (Call copy number aberrations from aCGH log ratios or NGS data.)
# INPUT segmented.tsv: segmented.tsv TYPE GENE_EXPRS
# OUTPUT aberrations.tsv: aberrations.tsv
# OUTPUT aberration-frequencies.pdf: aberration-frequencies.pdf
# PARAMETER number.of.copy.number.states: "number of copy number states" TYPE [3: 3, 4: 4] DEFAULT 3 (Three states means calling loss vs. normal vs. gain, four states calls amplifications separately, and five also homozygous deletions.)
# PARAMETER organism: organism TYPE [human: human, other: other] DEFAULT human (Organism.)
# PARAMETER genome.build: "genome build" TYPE [GRCh37: GRCh37, NCBI36: NCBI36, NCBI35: NCBI35, NCBI34: NCBI34] DEFAULT GRCh37 (The genome build to use. GRCh37 = hg19, NCBI36 = hg18, NCBI35 = hg17, NCBI34 = hg16. Not used unless organism is set to human.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2013-04-12

source(file.path(chipster.common.path, 'CGHcallPlus.R'))

file <- 'segmented.tsv'
dat <- read.table(file, header=TRUE, sep='\t', quote='', row.names=1, as.is=TRUE, check.names=FALSE)

if (length(grep("^segmented\\.", names(dat))) == 0)
  stop('CHIPSTER-NOTE: This tool needs segmented input, so please first run the tool Segment copy number data.')

dat$chromosome[dat$chromosome=='X'] <- '23'
dat$chromosome[dat$chromosome=='Y'] <- '24'
dat$chromosome[dat$chromosome=='MT'] <- '25'
dat$chromosome <- as.integer(dat$chromosome)

logratios <- as.matrix(dat[,grep("^chip\\.", names(dat))])
segmented <- as.matrix(dat[,grep("^segmented\\.", names(dat))])

logratios[logratios==-Inf] <- -10
segmented[segmented==-Inf] <- -10

cgh.seg <- new('cghSeg', assayData=assayDataNew(copynumber=logratios, segmented=segmented), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))

remove <- rowSums(is.na(logratios)) + rowSums(is.na(segmented))
cgh.seg <- cgh.seg[remove == 0,]
cgh.seg <- cgh.seg[chromosomes(cgh.seg) < 24,]

# identify different matrices (chip, flag, ...) present in the data
x <- colnames(dat)
suffix <- sub('^chip\\.', '', x[grep('^chip\\.', x)[1]])
matrices <- sub(suffix, '', x[grep(suffix, x)])
annotations <- 1:length(x)
for (m in matrices) {
  annotations <- setdiff(annotations, grep(m, x))
}

dat3 <- dat[featureNames(cgh.seg), annotations]

chips <- colnames(dat)[grep("^chip\\.", names(dat))]

# clean up to save memory
rm(file, dat, logratios, segmented, remove, x, suffix, matrices, annotations, m)
gc()

cgh.psn <- postsegnormalize(cgh.seg)

rm(cgh.seg)
gc()

cgh.cal <- CGHcall(cgh.psn, nclass=as.integer(number.of.copy.number.states), organism=organism, build=genome.build, ncpus=4)
cgh <- ExpandCGHcall(cgh.cal, cgh.psn)

rm(cgh.psn, cgh.cal)
gc()

#dat3 <- fData(cgh)
#colnames(dat3) <- c('chromosome', 'start', 'end')
#dat3 <- dat3[featureNames(cgh),]

#for (col in c('cytoband', 'symbol', 'description', 'cnvs'))
#  if (col %in% colnames(dat))
#    dat3[,col] <- dat[rownames(dat3), col]

calls <- calls(cgh)

if (-2 %in% calls)
  dat3$del.freq <- round(rowMeans(calls == -2), digits=3)
dat3$loss.freq <- round(rowMeans(calls == -1), digits=3)
dat3$gain.freq <- round(rowMeans(calls == 1), digits=3)
if (2 %in% calls)
  dat3$amp.freq <- round(rowMeans(calls == 2), digits=3)

cols <- colnames(dat3)

dat3 <- cbind(dat3, calls)
cols <- c(cols, sub('chip.', 'flag.', chips))
rm(calls)
gc()

dat3 <- cbind(dat3, copynumber(cgh))
cols <- c(cols, chips)
gc()

dat3 <- cbind(dat3, segmented(cgh))
cols <- c(cols, sub('chip.', 'segmented.', chips))
gc()

if (number.of.copy.number.states=='5') {
  dat3 <- cbind(dat3, probdloss(cgh))
  cols <- c(cols, sub('chip.', 'probdel.', chips))
  gc()
}

dat3 <- cbind(dat3, probloss(cgh))
cols <- c(cols, sub('chip.', 'probloss.', chips))
gc()

dat3 <- cbind(dat3, probnorm(cgh))
cols <- c(cols, sub('chip.', 'probnorm.', chips))
gc()

dat3 <- cbind(dat3, probgain(cgh))
cols <- c(cols, sub('chip.', 'probgain.', chips))
gc()

if (number.of.copy.number.states=='4') {
  dat3 <- cbind(dat3, probamp(cgh))
  cols <- c(cols, sub('chip.', 'probamp.', chips))
  gc()
}

colnames(dat3) <- cols

dat3$chromosome <- as.character(dat3$chromosome)
dat3$chromosome[dat3$chromosome=='23'] <- 'X'
dat3$chromosome[dat3$chromosome=='24'] <- 'Y'
dat3$chromosome[dat3$chromosome=='25'] <- 'MT'

options(scipen=10)
write.table(dat3, file='aberrations.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# bitmap(file='aberration-frequencies.png', width=image.width/72, height=image.height/72)
pdf(file='aberration-frequencies.pdf', paper='a4r', width=0, height=0)
frequencyPlot(cgh)
dev.off()

# EOF
