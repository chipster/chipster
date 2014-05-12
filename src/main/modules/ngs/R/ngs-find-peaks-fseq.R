# TOOL ngs-find-peaks-fseq.R: "Find broad peaks using F-seq" (This tool will search for statistically significantly enriched broad peaks, such as regions of open chromatin or transcription factor / histone binding sites, in sequencing data. Peaks are identified by smoothing read count data to a continuous sequence density estimation-based probabilites that are directly proportional to the probability of seeing a sequence read at that location.) 
# INPUT alignment.txt: "Read file" TYPE GENERIC 
# OUTPUT peaks.bed: "True enriched peaks in a format compatible with the Genome Browser"
# PARAMETER file.format: "File format" TYPE [ELAND, SAM, BAM, BED] DEFAULT BAM (The format of the input files.)
# PARAMETER score.threshold: "Score cutoff" TYPE DECIMAL FROM 0 TO 100 DEFAULT 4 (The cutoff for statistical significance. Larger values reduce the false discovery rate and produce shorter peaks lists.)
# PARAMETER fragment.size: "Fragment size" TYPE INTEGER FROM -1 TO 8000 DEFAULT -1 (Fragment size of the sequencing library. Set to -1, if your data contains more than 50000 mappable reads, in which case the fragment size is inferred from data.)
# PARAMETER feature.size: "Feature length" TYPE INTEGER FROM 0 TO 8000 DEFAULT 800 (The estimated feature length. The parameter controls the smoothness of the kernel density estimates. Larger values will lead to smoother kernel density estimation. In the case of DNase-seq, set this parameter to 0.)
# PARAMETER extend.reads: "Extend alignments" TYPE INTEGER FROM 0 TO 8000 DEFAULT 0 (Artificially extend each mapped read to a desired length, typically to the mean fragment length.)

# MK, 20.3.2012                                      
# F-seq is a tool to find broad peaks from pre-aligned SAM/BAM/ELAND files. The analysis needs to  #
# be done for each file separately. Found peaks are typically binned into categories representing  #
# exons (first, last all), introns, TSS, promoters, intergenic. Motifs can also be searched        #

#Export2sam.pl is a part of samtools and is used to convert ELAND outputs to SAM
eland_to_sam.binary <- c(file.path(chipster.tools.path, "samtools-0.1.18", "misc", "export2sam.pl"))
#samtools is used to convert SAM to BAM
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
#bedtools is used to convert BAM to BED
bam_to_bed.binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "bamToBed"))
#In-house perl script to extend reads from 5'
read_ex.binary <- c(file.path(chipster.tools.path, "fseq", "bin", "read_extend_bed.pm"))
#Binary file for F-seq which is used to call peaks
fseq.binary <- c(file.path(chipster.tools.path, "fseq", "bin", "fseq"))

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
	file.format <- "BAM"
	input.file <- "alignment.bam"
} 	

if (file.format == "BAM") {
	system(paste(bam_to_bed.binary, "-i ", input.file, ">alignment.bed"))
	input.file <- "alignment.bed"
} 

#In some articles, reads have been artifically extended to full fragment length. 
if (extend.reads > 0) {
	system(paste(read_ex.binary, " ",  input.file, " ", extend.reads, ">alignment.extend.bed"))
	input.file <- "alignment.extend.bed"
	#system(paste("cp -r ", getwd(), "/tmp/xxx/."))
}

#Count number of aligned reads from the BED-file
data <- read.table(file=input.file, sep="\t");
nread <- nrow(data);

#Output of F-seq is printed to a folder that must exists before F-seq is executed
dir.create("alignment_fseq")

#Execution of F-seq. Human data could benefit from background files that filter off uncertain areas. 
#However, no script exists for creating background files at the moment. Default background files may be suboptimal
if(fragment.size == -1) {
	if(nread < 50000) { stop("CHIPSTER-NOTE: F-seq needs at least 50000 aligned reads for fragment size estimation")  }
	system(paste(fseq.binary, "-l", feature.size," -t ", score.threshold," -of bed -o alignment_fseq", input.file));	
} else {
	system(paste(fseq.binary, "-f", fragment.size," -l", feature.size," -t ", score.threshold," -of bed -o alignment_fseq", input.file));	
}

#Output of F-seq is a bunch of files stored in folder -o. 
files 	<- list.files(path="alignment_fseq", pattern=".bed")
data 	<- NULL;
for(i in 1:length(files)) {
	chr_file <- files[i];
	chr_data <- read.table(paste("alignment_fseq/",chr_file, sep=""), sep="\t")
	data <- rbind(data, chr_data)
}
colnames(data) <- c("chr", "start", "end", "ID", "score")
rownames(data) <- data[,4]

#Sorting the BED
source(file.path(chipster.common.path, "bed-utils.R"))
if (nrow(data) > 1){	
	data <- sort.bed(data)
}
write.table(data, file="peaks.bed", sep="\t", row.names=F, col.names=F, quote=F)