# TOOL bwa-mem.R: "BWA MEM for single or paired end reads" (Aligns reads to genomes using the BWA MEM algorithm. If just one read file is given, then a single end analysis is run. If two read files are given, then mapping is done in paired end mode. Results are sorted and indexed BAM files, which are ready for viewing in the Chipster genome browser. 
# Note that this BWA tool uses publicly available genomes. If you would like to align reads against your own reference genome, please use the tool \"BWA MEM for single or paired end data with own genome\".)
# INPUT reads.fastq: "Read file 1" TYPE GENERIC 
# INPUT OPTIONAL reads2.fastq: "Read file 2 for paired end reads" TYPE GENERIC 
# OUTPUT bwa.bam 
# OUTPUT bwa.bam.bai 
# OUTPUT bwa.log 
# PARAMETER organism: "Organism" TYPE [Arabidopsis_thaliana.TAIR10.26, Bos_taurus.UMD3.1, Canis_familiaris.CanFam3.1, Drosophila_melanogaster.BDGP6, Felis_catus.Felis_catus_6.2, Gallus_gallus.Galgal4, Gasterosteus_aculeatus.BROADS1, Halorubrum_lacusprofundi_atcc_49239.GCA_000022205.1.26, Homo_sapiens.GRCh37.75, Homo_sapiens.GRCh38, Homo_sapiens_mirna, Medicago_truncatula.GCA_000219495.2.26, Mus_musculus.GRCm38, Mus_musculus_mirna, Ovis_aries.Oar_v3.1, Populus_trichocarpa.JGI2.0.26, Rattus_norvegicus_mirna, Rattus_norvegicus.Rnor_5.0, Schizosaccharomyces_pombe.ASM294v2.26, Sus_scrofa.Sscrofa10.2, Vitis_vinifera.IGGP_12x.26, Yersinia_enterocolitica_subsp_palearctica_y11.GCA_000253175.1.26] DEFAULT Homo_sapiens.GRCh38 (Genome that you would like to align your reads against.)
# PARAMETER OPTIONAL minseedlen: "Minimum seed length" TYPE INTEGER DEFAULT 19 (Matches shorter than this will be missed when looking for maximal exact matches or MEMs in the first alignment phase.)
# PARAMETER OPTIONAL bandwith: "Maximum gap length" TYPE INTEGER DEFAULT 100 (Gaps longer than this will not be found. Note also scoring matrix and hit length affect the maximum gap length, in addition to this band width parameter.)
# PARAMETER OPTIONAL matchscore: "Match score" TYPE INTEGER DEFAULT 1 (Score for a matching base.)
# PARAMETER OPTIONAL mismatchscore: "Mismatch penalty" TYPE INTEGER DEFAULT 4 (Penalty for a mismatching base.)
# PARAMETER OPTIONAL gapopen: "Gap opening penalty" TYPE INTEGER DEFAULT 6 (Gap opening penalty.)
# PARAMETER OPTIONAL gapextension: "Gap extension penalty" TYPE INTEGER DEFAULT 1 (Gap extension penalty.)
# PARAMETER OPTIONAL clippenalty: "Penalty for end clipping" TYPE INTEGER DEFAULT 5 (Penalty for 5\'- and 3\'-end clipping. When performing the Smith-Waterman extension of the seed alignments, BWA-MEM keeps track of the best score reaching the end of the read. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied.) 
# PARAMETER OPTIONAL rgid: "Read group dentifier" TYPE STRING (Read group identifier. If you want to add the read group line in the BAM file, you have to give this information.)
# PARAMETER OPTIONAL rgsm: "Sample name for read group" TYPE STRING (The name of the sample sequenced in this read group. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL rgpl: "Platform for read group" TYPE [ none: "Not defined", ILLUMINA, SOLID, LS454, HELICOS, PACBIO] DEFAULT none (Platform\/technology used to produce the read. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL rglb: "Library identifier for read group" TYPE STRING (DNA preparation library identifier. The Mark Duplicates tool uses this field to determine which read groups might contain molecular duplicates, in case the same DNA library was sequenced on multiple lanes. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)

# KM 28.08.2015

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")

pe_check.command <- ("ls reads2.fastq 2>/dev/null | wc -l")

paired.end <- system(pe_check.command, intern = TRUE )

if ( paired.end == 1 ){
	unzipIfGZipFile("reads2.fastq")
}

# bwa
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa mem")
bwa.genome <- file.path(chipster.tools.path, "genomes", "indexes", "bwa", organism)
#command.start <- paste("bash -c '", bwa.binary)
command.start <-(bwa.binary)

bwa.parameters <- paste("-M", "-k", minseedlen, "-w", bandwith, "-A", matchscore, "-B", mismatchscore, "-O", gapopen, "-E", gapextension, "-L",  clippenalty )

#Read group definition
if ( nchar(rgid) > 0 ){
	rg.string <-("'@RG")
	rg.string <- paste(rg.string, "\\tID:", rgid, sep="")
	if ( nchar(rgsm) > 0 ){
	   rg.string <- paste(rg.string, "\\tSM:", rgsm, sep="")
    }
	if ( rgpl != "none" ){
    	rg.string <- paste(rg.string, "\\tPL:", rgpl, sep="")
	}
	if ( nchar(rglb) > 0 ){
    	rg.string <- paste(rg.string, "\\tLB:", rglb, sep="")
	}
	rg.string <- paste(rg.string, "'", sep="")
	bwa.parameters <- paste(bwa.parameters,  "-R", rg.string )
}

# command ending
command.end <- paste(bwa.genome, "reads.fastq 1> alignment.sam 2>> bwa.log")
if ( paired.end == 1 ){
	command.end <- paste(bwa.genome, "reads.fastq reads2.fastq 1> alignment.sam 2>> bwa.log")
}	
	
bwa.command <- paste(command.start, bwa.parameters, command.end)

echo.command <- paste("echo '", bwa.binary, bwa.parameters, bwa.genome, "reads.fastq ' > bwa.log" )
if ( paired.end == 1 ){
  echo.command <- paste("echo '", bwa.binary, bwa.parameters, bwa.genome, "reads.fastq raeds2.fastq ' > bwa.log" )
}

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

# Handle output names
#
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

# Determine base name
basename <- strip_name(inputnames$reads.fastq)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("bwa.bam", paste(basename, ".bam", sep =""))
outputnames[2,] <- c("bwa.bam.bai", paste(basename, ".bam.bai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)

