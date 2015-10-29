#!/bin/bash
#
# K.M. 29.9. 2015
TOOLS_PATH=/opt/chipster/tools


#a5_miseq (päällekkäisiä työkaluja)
#PATH

#bedtools
BIN_PATH=bedtools/bin

export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#blast
BIN_PATH=blast/bin
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH



#bowtie
BIN_PATH=bowtie
export PATH=${PATH}::$TOOLS_PATH/$BIN_PATH

#bowtie2
BIN_PATH=bowtie2
export PATH=${PATH}::$TOOLS_PATH/$BIN_PATH

#bwa
BIN_PATH=bwa
export PATH=${PATH}::$TOOLS_PATH/$BIN_PATH

#ClusterBuster
BIN_PATH=ClusterBuster
export PATH=${PATH}::$TOOLS_PATH/$BIN_PATH

#cufflinks
BIN_PATH=cufflinks
export PATH=${PATH}::$TOOLS_PATH/$BIN_PATH

#cufflinks2
BIN_PATH=cufflinks2
export PATH=${PATH}::$TOOLS_PATH/$BIN_PATH

#dexseq-exoncounts
DEXSEQ_PATH=$TOOLS_PATH/dexseq-exoncounts
alias dexseq_count='python $DEXSEQ_PATH/dexseq_count.py' 

#dimont
DIMONT_PATH=$TOOLS_PATH/dimont
alias dimont='java -jar $DIMONT_PATH/dimond'

#edirect
BIN_PATH=edirect
export PATH=${PATH}::$TOOLS_PATH/$BIN_PATH

#emboss
BIN_PATH=emboss/bin
export PATH=${PATH}::$TOOLS_PATH/$BIN_PATH

#express
BIN_PATH=express
export PATH=${PATH}::$TOOLS_PATH/$BIN_PATH

#FastQC
BIN_PATH=FastQC
export PATH=${PATH}::$TOOLS_PATH/$BIN_PATH

#fastx
# on /usr/bin

#fseq
BIN_PATH=fseq/bin
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#fseq_bff

#GenomeAnalysisTK2
GATK2_PATH=$TOOLS_PATH/GenomeAnalysisTK2
alias gatk='java -jar $GATK2_PATH/GenomeAnalysisTKLite.jar'

#genomes
GENOMES_FASTA=$TOOLS_PATH/genomes/fasta
GENOMES_GTF=$TOOLS_PATH/genomes/fasta
GENOMES_BOWTIE=$TOOLS_PATH/genomes/indexes/bowtie
GENOMES_BOWTIE2=$TOOLS_PATH/genomes/indexes/bowtie2
GENOMES_BWA=$TOOLS_PATH/genomes/indexes/bwa
GENOMES_TOPHAT=$TOOLS_PATH/genomes/indexes/tophat

#htseq
# on /usr/local/bin

#macs
# on /usr/local/bin

#mafft
BIN_PATH=mafft/bin
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#meme
BIN_PATH=meme/bin
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#miRNA_mappings

#mothur
BIN_PATH=mothur
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#mothur-data

#picard-tools

#primer3
BIN_PATH=primer3/src/
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#prinseq
BIN_PATH=prinseq
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#QDNAseq

#R
BIN_PATH=R-3.1.2/bin
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#RSeQC
#samtools
BIN_PATH=samtools
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#seqtk
BIN_PATH=seqtk
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#sratoolkit
BIN_PATH=sratoolkit/bin
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#tabix
BIN_PATH=tabix
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#tagcleaner
TC_PATH=$TOOLS_PATH/tagcleaner
alias tagcleaner='perl $TC_PATH/tagcleaner.pl'

#tophat
BIN_PATH=tophat
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

#tophat2
BIN_PATH=tophat2
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH
#trimmomatic
TRIMMOMATIC_PATH=$TOOLS_PATH/trimmomatic
alias trimmomatic='java -jar $TRIMMOMATIC_PATH/trimmomatic-0.32.ja'

#vcftools
BIN_PATH=vcftools/bin
export PATH=${PATH}:$TOOLS_PATH/$BIN_PATH

