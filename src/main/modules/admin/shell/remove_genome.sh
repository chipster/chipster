#!/bin/bash

species=x
chipster_path=x 
while [[ $# -ge 1 ]]
do
  case "$1" in
              '-chipster_path')
	      chipster_path="$2"
                shift
                shift
              ;;
              #
              '-species')
                species="$2"
                shift 
                shift 
              ;;
              *)
                species="$1"
                shift 
              ;; 
  esac
done

if [[ $chipster_path == "x" ]] 
then
  echo "Please define the location of your chipster installation with option:"
  echo "-chipster_path"
  exit 1
fi


if [[ $species == "x" ]] 
then
  echo "Available species names:" 
  ls "$chipster_path"/tools/genomes/fasta/nochr/*.fa "$chipster_path"/tools/genomes/fasta/*.fa | awk -F "/" '{ print $NF}'   | sort | uniq | sed s/".toplevel.fa"/".toplevel"/g  
  exit 1
fi

tools_path="$chipster_path""/tools"
export PATH=${PATH}:"$chipster_path"/comp/modules/admin/shell/:"$tools_path"/emboss/bin/

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


#remove bwa_indexes
if [[ $INDEX_BWA -eq 1 ]]
then
  echo -------------------------------------------------------
  echo removing BWA indexes for $species
  rm $tools_path/bwa_indexes/$species.fa.amb
  rm $tools_path/bwa_indexes/$species.fa.ann
  rm $tools_path/bwa_indexes/$species.fa.bwt
  rm $tools_path/bwa_indexes/$species.fa.pac 
  rm $tools_path/bwa_indexes/$species.fa.sa
fi

#remove bowtie_indexes
if [[ $INDEX_BOWTIE -eq 1 ]]
then
  echo -------------------------------------------------------
  echo removing Bowtie indexes for $species
  rm $tools_path/bowtie/indexes/$species.1.ebwt
  rm $tools_path/bowtie/indexes/$species.2.ebwt
  rm $tools_path/bowtie/indexes/$species.3.ebwt
  rm $tools_path/bowtie/indexes/$species.4.ebwt
  rm $tools_path/bowtie/indexes/$species.rev.1.ebwt
  rm $tools_path/bowtie/indexes/$species.rev.2.ebwt
fi

#remove bowtie2_indexes
if [[ $INDEX_BOWTIE2 -eq 1 ]]
then
  echo -------------------------------------------------------
  echo removing Bowtie2 indexes for $species
  rm $tools_path/bowtie2/indexes/$species.1.bt2
  rm $tools_path/bowtie2/indexes/$species.2.bt2
  rm $tools_path/bowtie2/indexes/$species.3.bt2
  rm $tools_path/bowtie2/indexes/$species.4.bt2
  rm $tools_path/bowtie2/indexes/$species.rev.1.bt2
  rm $tools_path/bowtie2/indexes/$species.rev.2.bt2
fi


#remove fasta files
if [[ -e $tools_path/genomes/fasta/$species.fa ]]
then
   rm $tools_path/genomes/fasta/$species.fa*
fi 

if [[ -e $tools_path/genomes/fasta/nochr/$species.fa ]]
then 
   rm $tools_path/genomes/fasta/nochr/$species.fa*
fi 

#Remove gtf fiiles
species_gtf=$(echo $species | sed s/".dna.toplevel"/""/g) 
if [[ -e  $tools_path/genomes/gtf/$species_gtf.gtf ]]
then
  rm $tools_path/genomes/gtf/$species_gtf.*
fi

#Remove genomebrowser files.
version=$(echo $species | awk -F "." '{for (i=2; i<=NF; i++) printf $i"." }' | awk -F ".dna." '{print $1}')
species_name=$(echo $species | awk -F "." '{print $1}')
rm -rf $tools_path/genomebrowser/annotations/$species_name/$version



awk '{ if ( $3 != "'$species.fa'" ) print $0}' $tools_path/genomes/genome_list > $tools_path/genomes/genome_list.tmp.$$
rm $tools_path/genomes/genome_list
mv $tools_path/genomes/genome_list.tmp.$$ $tools_path/genomes/genome_list


