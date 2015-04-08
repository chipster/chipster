# TOOL bwa-mem-paired-end.R: "BWA MEM for paired-end reads" (BWA aligns reads to genomes and transcriptomes using BWA MEM algorithm. Results are sorted and indexed bam files, which are ready for viewing in the Chipster genome browser. 
# Note that this BWA tool uses publicly available genomes. If you would like to align reads against your own datasets, please use the tool \"BWA MEM for single end reads and own genome\".)
# INPUT reads1.txt: "Paired-end read set 1 to align" TYPE GENERIC 
# INPUT reads2.txt: "Paired-end read set 2 to align" TYPE GENERIC 
# OUTPUT bwa.bam 
# OUTPUT bwa.bam.bai 
# OUTPUT bwa.log 
# PARAMETER organism: "Genome or transcriptome" TYPE [Arabidopsis_thaliana.TAIR10.26, Bos_taurus.UMD3.1, Canis_familiaris.CanFam3.1, Drosophila_melanogaster.BDGP6, Felis_catus.Felis_catus_6.2, Gallus_gallus.Galgal4, Gasterosteus_aculeatus.BROADS1, Halorubrum_lacusprofundi_atcc_49239.GCA_000022205.1.26, Homo_sapiens.GRCh37.75, Homo_sapiens.GRCh38, Homo_sapiens_mirna, Medicago_truncatula.GCA_000219495.2.26, Mus_musculus.GRCm38, Mus_musculus_mirna, Ovis_aries.Oar_v3.1, Populus_trichocarpa.JGI2.0.26, Rattus_norvegicus_mirna, Rattus_norvegicus.Rnor_5.0, Schizosaccharomyces_pombe.ASM294v2.26, Sus_scrofa.Sscrofa10.2, Vitis_vinifera.IGGP_12x.26, Yersinia_enterocolitica_subsp_palearctica_y11.GCA_000253175.1.25] DEFAULT Homo_sapiens.GRCh38 (Genome or transcriptome that you would like to align your reads against.)
# PARAMETER mode: "Data source" TYPE [normal: " Illumina, 454, IonTorrent reads longer than 70 base pairs", pacbio: "PacBio subreads"] DEFAULT normal (Defining the type of reads will instruct the tool to use a predefined set of parameters optimized for that read type.)

# KM 11.11.2014

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads1.txt")
unzipIfGZipFile("reads2.txt")

# bwa
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa mem")
bwa.genome <- file.path(chipster.tools.path, "genomes", "indexes", "bwa", organism)
command.start <- paste("bash -c '", bwa.binary)

mode.parameters <- ifelse(mode == "pacbio", "-x pacbio", "")

# command ending
command.end <- paste(bwa.genome, "reads1.txt reads2.txt 1> alignment.sam 2>> bwa.log'")

# run bwa alignment
bwa.command <- paste(command.start, mode.parameters, command.end)

echo.command <- paste("echo '", bwa.binary , mode.parameters, bwa.genome, "reads1.txt reads2.txt' > bwa.log" )
#stop(paste('CHIPSTER-NOTE: ', bwa.command))
system(echo.command)
system(bwa.command)
		
# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# convert sam to bam
system(paste(samtools.binary, "view -bS alignment.sam -o alignment.bam"))

# sort bam
system(paste(samtools.binary, "sort alignment.bam alignment.sorted"))

# index bam
system(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
system("mv alignment.sorted.bam bwa.bam")
system("mv alignment.sorted.bam.bai bwa.bam.bai")
