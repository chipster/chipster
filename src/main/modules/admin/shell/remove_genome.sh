#!/bin/bash 

tools_path="/opt/chipster/tools/"
export PATH=${PATH}:/opt/chipster4/comp/modules/admin/shell/:/opt/chipster/tools/emboss/bin/

ensembl=0
fasta=0
gtf=0
text=0
length=0
version="0.0"
location=$(pwd)
INDEX_BWA=1
INDEX_BOWTIE=1
INDEX_BOWTIE2=1

species="$1"

#make bwa_indexes
if [[ $INDEX_BWA -eq 1 ]]
then
  echo -------------------------------------------------------
  echo removing BWA indexes for $genome_fasta
  rm $tools_path/bwa_indexes/$species.fa.amb
  rm $tools_path/bwa_indexes/$species.fa.ann
  rm $tools_path/bwa_indexes/$species.fa.bwt
  rm $tools_path/bwa_indexes/$species.fa.pac 
  rm $tools_path/bwa_indexes/$species.fa.sa
fi

#make bowtie_indexes
if [[ $INDEX_BOWTIE -eq 1 ]]
then
  echo -------------------------------------------------------
  echo removing Bowtie indexes for $genome_fasta
  rm $tools_path/bowtie/indexes/$species.1.ebwt
  rm $tools_path/bowtie/indexes/$species.2.ebwt
  rm $tools_path/bowtie/indexes/$species.3.ebwt
  rm $tools_path/bowtie/indexes/$species.4.ebwt
  rm $tools_path/bowtie/indexes/$species.rev.1.ebwt
  rm $tools_path/bowtie/indexes/$species.rev.2.ebwt
fi

#make bowtie_indexes
if [[ $INDEX_BOWTIE2 -eq 1 ]]
then
  echo -------------------------------------------------------
  echo removing Bowtie2 indexes for $genome_fasta
  rm $tools_path/bowtie2/indexes/$species.1.bt2
  rm $tools_path/bowtie2/indexes/$species.2.bt2
  rm $tools_path/bowtie2/indexes/$species.3.bt2
  rm $tools_path/bowtie2/indexes/$species.4.bt2
  rm $tools_path/bowtie2/indexes/$species.rev.1.bt2
  rm $tools_path/bowtie2/indexes/$species.rev.2.bt2
fi


rm $tools_path/genomes/fasta/$species.fa
rm $tools_path/genomes/fasta/nochr/$species.fa

species_gtf=$(echo $species | sed s/".dna.toplevel"/""/g) 
rm $tools_path/genomes/gtf/$species_gtf.gtf



awk '{ if ( $2 != "'$species'" ) print $0}' $tools_path/genomes/genome_list > $tools_path/genomes/genome_list.tmp.$$
rm $tools_path/genomes/genome_list
mv $tools_path/genomes/genome_list.tmp.$$ $tools_path/genomes/genome_list


