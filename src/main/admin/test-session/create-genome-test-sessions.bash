#!/bin/bash

# exit on error
set -e

bash run-tool-for-all-genomes.bash  h1-hESC_RNAseq.fastq bowtie.R bowtie-indexes.zip
bash run-tool-for-all-genomes.bash  h1-hESC_RNAseq.fastq bowtie2.R bowtie2-indexes.zip
bash run-tool-for-all-genomes.bash  h1-hESC_RNAseq.fastq tophat2.R tophat2-indexes.zip
bash run-tool-for-all-genomes.bash  h1-hESC_RNAseq.fastq bwa.R bwa-indexes.zip
#bash run-tool-for-all-genomes.bash  h1-hESC.bam DEXSeq.R DEXSeq-files.zip
bash run-tool-for-all-genomes.bash  h1-hESC.bam cufflinks2.R GTF-and-FASTA-files.zip

mv bowtie.R.genomes baseline
genome_diff=$(diff --from-file baseline *.genomes)
if [ "$genome_diff" ]
then
  echo "Warning: tools have different genome lists"
  echo ""
  echo "$genome_diff"
  exit 1
fi

rm *.genomes
