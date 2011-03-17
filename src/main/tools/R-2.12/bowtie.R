# TOOL "NGS" / bowtie.sadl: Bowtie (How to tie a bow.) 
# INPUT reads.txt: "Reads to align" TYPE GENERIC
# OUTPUT alignment-bowtie.bam
# OUTPUT alignment-bowtie.bam.bai
# OUTPUT bowtie.log
# PARAMETER genome: "Genome" TYPE [m_musculus_ncbi37: "m_musculus_ncbi37", hg19: "hg19", mm9: "mm9", rn4: "rn4"] DEFAULT m_musculus_ncbi37 (Genome.)
# PARAMETER max.mismatches: "Max mismatches in seed" TYPE [0, 1, 2, 3] DEFAULT 2 (Max mismatches in seed.)

# run bowtie
bowtie.binary <- c(file.path(chipster.tools.path, "bowtie", "bowtie"))
bowtie.command <- paste("bash -c '", bowtie.binary, "-n", max.mismatches, "-q --best -S", genome, "reads.txt 1> alignment.sam 2> bowtie.log'")
system(bowtie.command)

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# convert sam to bam
system(paste(samtools.binary, "view -bS alignment.sam -o alignment.bam"))

# sort bam
system(paste(samtools.binary, "sort alignment.bam alignment.sorted"))

# index bam
system(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
system("mv alignment.sorted.bam alignment-bowtie.bam")
system("mv alignment.sorted.bam.bai alignment-bowtie.bam.bai")