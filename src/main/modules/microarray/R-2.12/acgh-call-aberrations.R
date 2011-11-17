# TOOL acgh-call-aberrations.R: "Call copy number aberrations from aCGH data" (Call copy number aberrations from aCGH log ratios.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT aberrations.tsv: aberrations.tsv 
# OUTPUT aberrations.png: aberrations.png 
# PARAMETER normalization: normalization TYPE [median: median, mode: mode, none: none] DEFAULT none (Normalization method.)
# PARAMETER number.of.chromosomes: number.of.chromosomes TYPE INTEGER DEFAULT 23 (Number of chromosomes. Usually 23 for sex-matched reference samples and 22 otherwise.)
# PARAMETER number.of.copy.number.states: number.of.copy.number.states TYPE [3: 3, 4: 4] DEFAULT 3 (Whether to call loss vs. normal vs. gain or loss vs. normal vs. gain vs. amplification.)
# PARAMETER minimum.number.of.probes.per.segment: minimum.number.of.probes.per.segment TYPE [2: 2, 3: 3, 4: 4, 5: 5] DEFAULT 2 (Minimum number of probes per segment.)
# PARAMETER minimum.number.of.sds.between.segments: minimum.number.of.sds.between.segments TYPE DECIMAL FROM 0 TO 10 DEFAULT 0 (Minimum number of standard deviations required between segments.)
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image.)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image.)
# PARAMETER organism: organism TYPE [human: human, other: other] DEFAULT human (Organism.)
# PARAMETER genome.build: genome.build TYPE [GRCh37: GRCh37, NCBI36: NCBI36, NCBI35: NCBI35, NCBI34: NCBI34] DEFAULT GRCh37 (The genome build to use. GRCh37 = hg19, NCBI36 = hg18, NCBI35 = hg17, NCBI34 = hg16. Not used unless organism is set to human.)

# detect-copy-number-aberrations.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2011-11-17

library(CGHcall) # source CGHcallPlus

dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

pos <- c('chromosome','start','end')
if (length(setdiff(pos, colnames(dat)))!=0)
  stop('CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.')

dat2 <- data.frame(probe=rownames(dat), dat[,c('chromosome', 'start', 'end')], dat[,grep("chip", names(dat))], stringsAsFactors=FALSE)
chips <- colnames(dat)[grep("chip", names(dat))]

if (ncol(dat2)==4)
  stop('CHIPSTER-NOTE: No array data found. The input file must have columns labeled with "chip.".')
if (ncol(dat2)==5)
  colnames(dat2)[5] <- chips[1]

dat2$chromosome[dat2$chromosome=='X'] <- 23
dat2$chromosome[dat2$chromosome=='Y'] <- 24
dat2$chromosome[dat2$chromosome=='MT'] <- 25
dat2$chromosome <- as.integer(dat2$chromosome)

cgh.raw <- make_cghRaw(dat2)
cgh.pre <- preprocess(cgh.raw, nchrom=number.of.chromosomes)
cgh.nor <- normalize(cgh.pre, method=normalization)
cgh.seg <- segmentData(cgh.nor, min.width=as.integer(minimum.number.of.probes.per.segment), undo.splits='sdundo', undo.SD=minimum.number.of.sds.between.segments)
cgh.psn <- postsegnormalize(cgh.seg)
# cgh.cal <- CGHcall(cgh.psn, nclass=as.integer(number.of.copy.number.states), organism=organism, build=genome.build)
cgh.cal <- CGHcall(cgh.psn, nclass=as.integer(number.of.copy.number.states), organism=organism) # change to CGHcallPlus
cgh <- ExpandCGHcall(cgh.cal, cgh.psn)

dat3 <- data.frame(cgh@featureData@data)
colnames(dat3) <- c('chromosome', 'start', 'end')

for (col in c('cytoband', 'symbol', 'description', 'cnvs'))
  if (col %in% colnames(dat))
    dat3[,col] <- dat[rownames(dat3), col]

calls <- assayDataElement(cgh, 'calls')

dat3$loss.freq <- round(rowMeans(calls == -1), digits=3)
dat3$gain.freq <- round(rowMeans(calls == 1), digits=3)
if (number.of.copy.number.states=='4' && 2 %in% assayDataElement(cgh, 'calls'))
  dat3$amp.freq <- round(rowMeans(calls == 2), digits=3)

colnames(calls) <- sub('chip.', 'flag.', chips)
dat3 <- cbind(dat3, calls)

copynumber <- assayDataElement(cgh, 'copynumber')
colnames(copynumber) <- chips
dat3 <- cbind(dat3, copynumber)

segmented <- assayDataElement(cgh, 'segmented')
colnames(segmented) <- sub('chip.', 'segmented.', chips)
dat3 <- cbind(dat3, segmented)

probloss <- assayDataElement(cgh, 'probloss')
colnames(probloss) <- sub('chip.', 'probloss.', chips)
dat3 <- cbind(dat3, probloss)

probnorm <- assayDataElement(cgh, 'probnorm')
colnames(probnorm) <- sub('chip.', 'probnorm.', chips)
dat3 <- cbind(dat3, probnorm)

probgain <- assayDataElement(cgh, 'probgain')
colnames(probgain) <- sub('chip.', 'probgain.', chips)
dat3 <- cbind(dat3, probgain)

if (number.of.copy.number.states=='4') {
  probamp <- assayDataElement(cgh, 'probamp')
  colnames(probamp) <- sub('chip.', 'probamp.', chips)
  dat3 <- cbind(dat3, probamp)
}

dat3$chromosome <- as.character(dat3$chromosome)
dat3$chromosome[dat3$chromosome=='23'] <- 'X'
dat3$chromosome[dat3$chromosome=='24'] <- 'Y'
dat3$chromosome[dat3$chromosome=='25'] <- 'MT'

write.table(dat3, file='aberrations.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

bitmap(file='aberrations.png', width=image.width/72, height=image.height/72)
# pdf(file='aberrations.pdf')
plot.summary(cgh) # change for CGHcallPlus
# frequencyPlot(cgh)
dev.off()

# EOF
