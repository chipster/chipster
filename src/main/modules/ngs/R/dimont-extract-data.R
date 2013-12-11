# TOOL dimontdata.R: "Dimont data extractor" (Extracts genomic regions specified in a BED-like file format in the annotated FastA format as required by Dimont.)
# INPUT regions.bed: "Genomic regions" TYPE GENERIC (The genomic regions to be extracted in a BED-like file format, e.g., BED, GTF, narrowPeak.)
# OUTPUT extracted.fa: "Extracted sequences" (The sequences extracted from the given genome using the supplied region specifications.)
# PARAMETER genome: "Genome" TYPE [hg19: "Human (hg19\)", mm10: "Mouse (mm10\)", rn4: "Rat (rn4\)", Halorubrum_lacusprofundi_ATCC_49239: "Halorubrum lacusprofundi ATCC 49239", Sus_scrofa.Sscrofa10.2.69.dna.toplevel: "Pig (sus_scrofa10.2.69\)",  athaliana.TAIR10: "A. thaliana (TAIR10\)", Arabidopsis_lyrata.v.1.0.16: "A. lyrata (1.0.16\)", ovis_aries_texel: "Sheep (oar3.1\)"] DEFAULT hg19 (Genome or transcriptome that the genomic regions refer to.)
# PARAMETER chromcol: "Chromosome column" TYPE INTEGER FROM 1 DEFAULT 1 (The column of the Regions file, which contains the chromosome information. Please note that row-names are not considered as a column.)
# PARAMETER startcol: "Start column" TYPE INTEGER FROM 1 DEFAULT 2 (The column of the Regions file containing the start position of the genomic region. Please note that row-names are not considered as a column.)
# PARAMETER seccol: "Meaning of second coordinate" TYPE [center: "Center of peak (relative to start\)", end: "End of peak (global coordinates\)"] DEFAULT end (The meaning of the second genomic coordinate. This may either be the position of the peak summit relative to the position in Start, or the end position of the peak.)
# PARAMETER seccoord: "Second coordinate" TYPE INTEGER FROM 1 DEFAULT 3 (The second genomic coordinate with meaning specified by parameter \"Meaning of second coordinate\". Please note that row-names are not considered as a column.)
# PARAMETER width: Width TYPE INTEGER FROM 1 DEFAULT 1000 (The width of the genomic region to be extracted. Recommended values: 1000 for ChIP-seq and 100 for ChIP-exo.)
# PARAMETER statcol: "Statistics column" TYPE INTEGER FROM 1 DEFAULT 5 (The column containing the peak statistics information (or another measure of peak confidence\).)

#MK: 11.12.2013. Since also row-names are printed, column indexes had to vbe shifted by one 

genome.fa<-file.path(chipster.tools.path, "genomes", "fasta", paste(genome,".fa",sep="",collapse=""))

tool<-file.path(chipster.tools.path,"dimont","extract_data_single_chipster.pl");

command<-paste("perl",tool,genome.fa,"regions.bed",(chromcol+1),(startcol+1),seccol,(seccoord+1),width,(statcol+1),"extracted.fa");

system(command)
