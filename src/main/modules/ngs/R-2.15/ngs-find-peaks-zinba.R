# TOOL ngs-find-peaks-zinba.R: "Find broad peaks using Zinba" (This tool provide means to find broad peaks, such as regions of open chromatin or transcription factor / histone binding sites from DNase, FAIRE or Chip-seq data. Please note that Zinba is extremely slow and analyses can last several hours.)
# INPUT alignment.txt: "Read file" TYPE GENERIC 
# OUTPUT peaks.bed: "True enriched peaks in a format compatible with the Genome Browser"
# PARAMETER file.format: "File format" TYPE [ELAND, SAM, BAM, BED] DEFAULT BAM (The format of the input files.)
# PARAMETER species: "Genome" TYPE [human, yeast] DEFAULT human (The species of the samples.)
# PARAMETER read.len: "Read length" TYPE [36, 50] DEFAULT 36 (The size of the read.)
# PARAMETER fragment.size: "Fragment length" TYPE INTEGER FROM 1 TO 8000 DEFAULT 300 (The estimated fragment length of the sequencing library.)
# PARAMETER broad.peaks: "Peak type" TYPE [broad, narrow] DEFAULT narrow (Peak type that is searched. If broad is selected, adjacent peaks are merged into a single peak over a much longer region. Otherwise only overlapping peaks are merged)
# PARAMETER max.hit: "Number of hits per read at maximum" TYPE INTEGER FROM 1 TO 100 DEFAULT 4 (Number of hits per read allowed during mapping process.)
# PARAMETER score.threshold: "FDR cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (The FDR cutoff for statistical significance. Smaller values reduce the false discovery rate and produce shorter peaks lists.)
# PARAMETER model.chr: "One chromosome to use for model selection" TYPE STRING DEFAULT chr22 (One chromosome name to use for model selection.)


######################################################
#                                                    #
# MK, 20.3.2012                                      #
#                                                    #
# Zinba is a tool to find broad peaks from pre-      #
# aligned SAM/BAM/ELAND files. The analysis needs to #
# be done for each file separately. Peaks are        #
# typically binned into categories representing      #
# exons (first, last all), introns, TSS, promoters,  #
# intergenic. Motifs can be searched from them also  #
#                                                    #
######################################################


#Export2sam.pl is a part of samtools and is used to convert ELAND outputs to SAM
eland_to_sam.binary <- c(file.path(chipster.tools.path, "samtools-0.1.18", "misc", "export2sam.pl"))
#samtools is used to convert SAM to BAM
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
#bedtools is used to convert BAM to BED
bam_to_bed.binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "bamToBed"))

#File conversions ending to BED
input.file <- "alignment.txt";
if(file.format == "ELAND") {
	system(paste(eland_to_sam.binary, "--read1=", input.file," > alignment.sam"))
	file.format <- "SAM"
	input.file  <- "alignment.sam"
}
if (file.format == "SAM") {
	system(paste(samtools.binary, "view -bS", input.file," -o alignmentUnsorted.bam"))
	system(paste(samtools.binary, "sort alignmentUnsorted.bam alignment"))
	system(paste(samtools.binary, "index alignment.bam"))
	file.format <- "BAM"
	input.file <- "alignment.bam"
} 	
if (file.format == "BAM") {
	system(paste(bam_to_bed.binary, "-i ", input.file, ">alignment.bed"))
	input.file <- "alignment.bed"
} 

#Load ZINBA
library(zinba)

#Modify arguments to fit Zinba
if(species=="human") 		{ species_nick  = "hg19" }
if(species=="mouse") 		{ species_nick  = "mm9" }
if(species=="arabidopsis") 	{ species_nick  = "tair10" }
if(species=="yeast") 		{ species_nick  = "sacCer3" }

if(broad.peaks=="broad") 	{ 
	winGap <- 5000
	selectcovs <- c("gcPerc", "align_perc")	
} else { 
	winGap <- 0
	selectcovs <- c("gcPerc", "exp_cnvwin_log", "align_perc")
}

#Zinba writes some files to the alignability folder. To fix this, new folder is created and softs link to 
#the real index files created
dir.name	<- paste("map", read.len,"_", species_nick, sep="");
index.dir	<- file.path(chipster.tools.path, "zinba", "map", dir.name)
dir.create(dir.name);
system(paste("ln -s", file.path(index.dir,"*"), " ", file.path(dir.name,".")));

#Generate your alignability directory
generateAlignability(
		mapdir=dir.name, 	    																						# mappability directory from unpacked mappability files
		twoBitFile=c(file.path(chipster.tools.path, "zinba", "2bit", paste(species_nick, ".2bit", sep=""))),            # path to downloaded genome build file in .2bit format
		outdir="zinba_background/",                                														# directory for processed files, used later in analysis
		athresh=max.hit,                                              													# number of hits per read allowed during mapping process
		extension=fragment.size,                                       													# average fragment library length
);

#Run the basealigncount function to generate the basecount file needed to obtain exact peak boundaries through 
#peak refinement.
basealigncount(
		inputfile=input.file, 																							# mapped sample reads
		outputfile="zinba_background/zinba_basecount.txt", 																						# output path
		extension=fragment.size, 																						# average fragment library length
		filetype="bed", 																								# either "bed", "bowtie", or "tagAlign"
		twoBitFile=c(file.path(chipster.tools.path, "zinba", "2bit", paste(species_nick, ".2bit", sep=""))), 			# path to downloaded genome build file in .2bit format
);

func <- function(this, ...) { 
	k<- list(...); 
	print(k); 
}

run.zinba(
		refinepeaks=1,          	              															# refine peaks? 1 for yes, 0 for no
		seq=input.file,   		                       														# path to mapped experimental reads
		input="none",               	               														# path to mapped input reads if available (default is "none")
		filetype="bed",	                                   													# either 'bed', 'bowtie', or 'tagAlign'
		threshold=score.threshold,                     														# FDR threshold, default is 0.05
		align="zinba_background/",                            												# path to alignability directory
		numProc=1,																							# number of CPUs to use, must be less than max available   (default 1)
		twoBit=c(file.path(chipster.tools.path, "zinba", "2bit", paste(species_nick, ".2bit", sep=""))),    # path to genome build in .2bit format
		outfile="zinba.out",  																				# prefix for outputted files
		extension=fragment.size,																			# average fragment library length (size selected)
		basecountfile="zinba_background/zinba_basecount.txt",												# path to basecount file if refinepeaks is 1
		printFullOut=0,   																					# print original data with enrichment estimates, 1 for yes (more space required), 0 for no (default)
		interaction=TRUE,   																				# whether or not to considering interaction during model selection, TRUE (default) or FALSE
		FDR=TRUE, 																							# either TRUE (default) or FALSE. If false, then uses posterior probability to threshold peaks using 1-threshold
		winGap=winGap,	
		selectmodel=T,
		selectchr=model.chr,		
		selectcovs=selectcovs,																				# vector of covariate names (characters) to consider in model selection  (can be  "gcPerc", "align_perc", "exp_cnvwin_log", and "input_count")
		pWinSize=200, 																						# sliding window size for local maximum detection (default 200 bp)
		winSize=250,
		offset=125,	
		cnvWinSize=100000,
		cnvOffset=2500,
		pquant=1,
		initmethod="count",
		diff=0,
		method="mixture",
		selecttype="dirty",
);


#system(paste("cp -r ", "/tmp/matti/ce671df2-8b95-4826-8259-950c658c2d77/zinba.out3.peaks.bed ", getwd(), "/zinba.out.peaks.bed", sep=""))
#system(paste("cp -r ", getwd(), "/tmp/matti/."))
#print(list.files())

#Read Output
data <- read.table(file="zinba.out.peaks.bed", sep="\t")
colnames(data) <- c("chr", "start", "end", "ID", "score", "strand")
rownames(data) <- paste(data[,4], 1:nrow(data), sep="_");

#Sorting the BED
source(file.path(gsub("R-2.15", "R-2.12", chipster.common.path), "bed-utils.R"))
if (nrow(data) > 1){	
	data <- sort.bed(data)
}
write.table(data, file="peaks.bed", sep="\t", row.names=F, col.names=F, quote=F)

