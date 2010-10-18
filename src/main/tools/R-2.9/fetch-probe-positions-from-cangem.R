# ANALYSIS "aCGH tools"/"Fetch probe positions from CanGEM" (Fetches microarray probe positions from the CanGEM database.)
# INPUT GENE_EXPRS normalized.tsv
# OUTPUT probe-positions.tsv
# PARAMETER platform [other, Affymetrix_GeneChip_Human_Mapping_50K_Xba, Affymetrix_Human_Genome_U133_Plus_2.0, Affymetrix_Human_Genome_U133A, Affymetrix_Human_Genome_U133A_2.0, Affymetrix_Human_Genome_U133B, Affymetrix_Human_Genome_U95A, Affymetrix_Human_Genome_U95Av2, Affymetrix_Human_Genome_U95B, Affymetrix_Human_Genome_U95C, Affymetrix_Human_Genome_U95D, Affymetrix_Human_Genome_U95E, Agilent_Human_1_cDNA_Microarray_G4100A, Agilent_Human_1A_Microarray_G4110A, Agilent_Human_1A_Microarray_V2_G4110B, Agilent_Human_1B_Microarray_G4111A, Agilent_Human_Genome_CGH_Oligo_Microarray_4x44K, Agilent_Human_Genome_CGH_Oligo_Microarray_Kit_185, Agilent_Human_Genome_CGH_Oligo_Microarray_Kit_244A_G4411B, Agilent_Human_Genome_CGH_Oligo_Microarray_Kit_44A_G4410A, Agilent_Human_Genome_CGH_Oligo_Microarray_Kit_44B_G4410B, Agilent_Human_miRNA_Microarray, Agilent_Human_miRNA_Microarray_V2, Agilent_Human_miRNA_Microarray_V3, Agilent_Human_Promoter_ChIP_on_Chip_Set_244K_Microarray_1_of_2, Agilent_Human_Promoter_ChIP_on_Chip_Set_244K_Microarray_2_of_2, Agilent_Oxford_2x105, Agilent_Whole_Human_Genome_Microarray_Kit_G4112A, Agilent_Whole_Human_Genome_Oligo_Microarray_4x44K, Turku_HUM_16K_cDNA, VUMC_Human_30K_60mer_oligo_array, Ensembl_Genes, Ensembl_Cytobands, Ensembl_Chromosomes] DEFAULT other (The microarray platform. If the one used is not on the list, but can be found from CanGEM, please select other and specify the accession number using the next parameter.)
# PARAMETER other.platform.accession STRING DEFAULT CG-PLM- (The accession of the platform. This is used only if other is selected in the previous parameter.)
# PARAMETER genome.build [GRCh37, NCBI36, NCBI35, NCBI34] DEFAULT GRCh37 (The genome build to use for adding the chromosome names and start and end base pair positions for the probes.)

# fetch-probe-positions-from-cangem.R
# Ilari Scheinin <firstname.lastname@gmail.com>
# 2010-10-12

# determine platform accession if a platform was chosen from the popup
if (platform != 'other') {
 cangem.platforms <- c(Affymetrix_GeneChip_Human_Mapping_50K_Xba='CG-PLM-18', Affymetrix_Human_Genome_U133_Plus_2.0='CG-PLM-10', Affymetrix_Human_Genome_U133A='CG-PLM-7', Affymetrix_Human_Genome_U133A_2.0='CG-PLM-9', Affymetrix_Human_Genome_U133B='CG-PLM-8', Affymetrix_Human_Genome_U95A='CG-PLM-11', Affymetrix_Human_Genome_U95Av2='CG-PLM-16', Affymetrix_Human_Genome_U95B='CG-PLM-12', Affymetrix_Human_Genome_U95C='CG-PLM-13', Affymetrix_Human_Genome_U95D='CG-PLM-14', Affymetrix_Human_Genome_U95E='CG-PLM-15', Agilent_Human_1_cDNA_Microarray_G4100A='CG-PLM-1', Agilent_Human_1A_Microarray_G4110A='CG-PLM-19', Agilent_Human_1A_Microarray_V2_G4110B='CG-PLM-20', Agilent_Human_1B_Microarray_G4111A='CG-PLM-24', Agilent_Human_Genome_CGH_Oligo_Microarray_4x44K='CG-PLM-27', Agilent_Human_Genome_CGH_Oligo_Microarray_Kit_185='CG-PLM-23', Agilent_Human_Genome_CGH_Oligo_Microarray_Kit_244A_G4411B='CG-PLM-17', Agilent_Human_Genome_CGH_Oligo_Microarray_Kit_44A_G4410A='CG-PLM-2', Agilent_Human_Genome_CGH_Oligo_Microarray_Kit_44B_G4410B='CG-PLM-4', Agilent_Human_miRNA_Microarray='CG-PLM-36', Agilent_Human_miRNA_Microarray_V2='CG-PLM-33', Agilent_Human_miRNA_Microarray_V3='CG-PLM-32', Agilent_Human_Promoter_ChIP_on_Chip_Set_244K_Microarray_1_of_2='CG-PLM-34', Agilent_Human_Promoter_ChIP_on_Chip_Set_244K_Microarray_2_of_2='CG-PLM-35', Agilent_Oxford_2x105='CG-PLM-31', Agilent_Whole_Human_Genome_Microarray_Kit_G4112A='CG-PLM-5', Agilent_Whole_Human_Genome_Oligo_Microarray_4x44K='CG-PLM-28', Turku_HUM_16K_cDNA='CG-PLM-3', VUMC_Human_30K_60mer_oligo_array='CG-PLM-30', Ensembl_Genes='CG-PLM-26', Ensembl_Cytobands='CG-PLM-6', Ensembl_Chromosomes='CG-PLM-37')
  other.platform.accession <- cangem.platforms[platform]
}

# check for valid platform accession
platform.accession <- toupper(other.platform.accession)
if (length(grep('^CG-PLM-[0-9]+$', platform.accession))==0)
  stop('CHIPSTER-NOTE: Not a valid platform accession: ', platform.accession)

dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE, row.names=1)

# remove probe positions if already present
dat$chromosome <- NULL
dat$start <- NULL
dat$end <- NULL
dat$cytoband <- NULL

# load platform
platform <- read.table(paste('http://www.cangem.org/download.php?platform=', platform.accession, '&flag=', genome.build, sep=''), sep='\t', header=TRUE, as.is=TRUE)
colnames(platform) <- tolower(colnames(platform))
colnames(platform)[colnames(platform)=='chr'] <- 'chromosome'
rownames(platform) <- platform[,1]
platform$chromosome <- factor(platform$chromosome, levels=c(1:22, "X", "Y", "MT"), ordered=TRUE)

dat2 <- cbind(platform[rownames(dat), c('chromosome', 'start', 'end')], dat, row.names=rownames(dat))
dat2 <- dat2[order(dat2$chromosome, dat2$start),]

write.table(dat2, file='probe-positions.tsv', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

# EOF