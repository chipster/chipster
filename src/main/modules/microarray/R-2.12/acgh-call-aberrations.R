# TOOL acgh-call-aberrations.R: "Call aberrations from segmented copy number data" (Call copy number aberrations from aCGH log ratios or NGS data.)
# INPUT segmented.tsv: segmented.tsv TYPE GENE_EXPRS
# OUTPUT aberrations.tsv: aberrations.tsv
# OUTPUT aberration-frequencies.pdf: aberration-frequencies.pdf
# PARAMETER number.of.copy.number.states: number.of.copy.number.states TYPE [3: 3, 4: 4] DEFAULT 3 (Whether to call loss vs. normal vs. gain or loss vs. normal vs. gain vs. amplification.)
# PARAMETER organism: organism TYPE [human: human, other: other] DEFAULT human (Organism.)
# PARAMETER genome.build: genome.build TYPE [GRCh37: GRCh37, NCBI36: NCBI36, NCBI35: NCBI35, NCBI34: NCBI34] DEFAULT GRCh37 (The genome build to use. GRCh37 = hg19, NCBI36 = hg18, NCBI35 = hg17, NCBI34 = hg16. Not used unless organism is set to human.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2012-01-05

source(file.path(chipster.tools.path, 'MPScall', 'CGHcallPlus-R-2.12.R'))

dat <- read.table('segmented.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

if (length(grep("^segmented\\.", names(dat))) == 0)
  stop('CHIPSTER-NOTE: This tool needs segmented input, so please first run the tool Segment copy number data.')

dat$chromosome[dat$chromosome=='X'] <- '23'
dat$chromosome[dat$chromosome=='Y'] <- '24'
dat$chromosome[dat$chromosome=='MT'] <- '25'
dat$chromosome <- as.integer(dat$chromosome)

logratios <- as.matrix(dat[,grep("^chip\\.", names(dat))])
segmented <- as.matrix(dat[,grep("^segmented\\.", names(dat))])

cgh.seg <- new('cghSeg', assayData=assayDataNew(copynumber=logratios, segmented=segmented), featureData=new('AnnotatedDataFrame', data=data.frame(Chromosome=dat$chromosome, Start=dat$start, End=dat$end, row.names=row.names(dat))))

chips <- colnames(dat)[grep("^chip\\.", names(dat))]

cgh.psn <- postsegnormalize(cgh.seg)
cgh.cal <- CGHcall(cgh.psn, nclass=as.integer(number.of.copy.number.states), organism=organism, build=genome.build)
cgh <- ExpandCGHcall(cgh.cal, cgh.psn)

dat3 <- data.frame(cgh@featureData@data)
colnames(dat3) <- c('chromosome', 'start', 'end')

for (col in c('cytoband', 'symbol', 'description', 'cnvs'))
  if (col %in% colnames(dat))
    dat3[,col] <- dat[rownames(dat3), col]

calls <- assayDataElement(cgh, 'calls')

dat3$loss.freq <- round(rowMeans(calls == -1), digits=3)
dat3$gain.freq <- round(rowMeans(calls == 1), digits=3)
if (number.of.copy.number.states=='4' && 2 %in% calls)
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

# bitmap(file='aberration-frequencies.png', width=image.width/72, height=image.height/72)
pdf(file='aberration-frequencies.pdf', paper='a4r', width=0, height=0)
frequencyPlot(cgh)
dev.off()

# EOF
