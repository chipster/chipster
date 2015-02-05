#!/bin/bash

if [[ $# -eq 0 ]]
then
  echo "Usage:"
  echo ""
  echo 'update_mirna "species name"'
  echo "" 
  exit 1
fi

species=$1
index_path=$2

INDEX_BOWTIE=(1)
INDEX_BOWTIE2=(1)
INDEX_BWA=(1)

#species="Homo sapiens"
species_name=$(echo $species | tr " " "_")


chipster_path=(/opt/chipster)
tools_path=($chipster_path/tools)
#index_path=($tools_path/genomes/indexes)

cd $tools_path/genomes/fasta

wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
gunzip mature.fa.gz

echo "Extracting mirnas for $species"
$chipster_path/comp/modules/misc/shell/textsearch_fasta.sh -emboss_path $tools_path/emboss/bin -sequence mature.fa -pattern "${species}" -outfile  ${species_name}_mirna_fasta_u

echo "------------------------------------------------------------------"
echo "results for mirna collection"
$chipster_path/comp/modules/misc/shell/infoseq_summary.sh $tools_path/emboss/bin ${species_name}_mirna_fasta_u
echo "------------------------------------------------------------------"

echo "Replacing U with T"
$tools_path/emboss/bin/biosed -targetregion U -replace T -sequence ${species_name}_mirna_fasta_u -outseq ${species_name}_mirna_fasta
 
genome=${species_name}_mirna

#make bwa_indexes
if [[ $INDEX_BWA -eq 1 ]]
then

  if [[ ! -e ${index_path}/bwa ]]
  then
    mkdir ${index_path}/bwa
  fi

  echo -------------------------------------------------------
  echo Calculating BWA indexes for ${genome}_fasta
  cd $index_path/bwa
  ln -s ../../fasta/${genome}_fasta ${genome}_fasta
  $tools_path/bwa/bwa index -p $genome  ${genome}_fasta
else
  echo "Skipping BWA indexing"
fi


#make bowtie_indexes
if [[ $INDEX_BOWTIE -eq 1 ]]
then
  if [[ ! -e ${index_path}/bowtie ]]
  then
    mkdir ${index_path}/bowtie
  fi

  echo -------------------------------------------------------
  echo Calculating Bowtie indexes for ${genome}_fasta 
  cd $index_path/bowtie
  ln -s ../../fasta/${genome}_fasta ${genome}_fasta
  $tools_path/bowtie/bowtie-build ${genome}_fasta $genome

  #check bowtie indexes
  n_genome=$(grep -c "^>" ${genome}_fasta)
  n_index=$($tools_path/bowtie/bowtie-inspect -n $genome | wc -l)

  if [[ $n_genome -ne $n_index ]]
  then
     echo "ERROR: Bowtie indexing of genome $genome failed"
     exit 1
  else
    echo Bowtie index of genome $genome OK
  fi
else
    echo "Skipping Bowtie indexing"
fi




#make bowtie2_indexes
if [[ $INDEX_BOWTIE2 -eq 1 ]]
then
  if [[ ! -e ${index_path}/bowtie2 ]]
  then
    mkdir ${index_path}/bowtie2
  fi


  echo -------------------------------------------------------
  echo "Calculating Bowtie2 indexes for ${genome}_fasta" 

  cd $index_path/bowtie2
  ln -s ../../fasta/${genome}_fasta ${genome}_fasta  
  ln -s ${genome}_fasta ${genome}.fa
  
  $tools_path/bowtie2/bowtie2-build ${genome}_fasta $genome

  #check bowtie2 indexes
  n_genome=$(grep -c "^>" ${genome}_fasta)
  n_index=$($tools_path/bowtie2/bowtie2-inspect -n $genome | wc -l)

  if [[ $n_genome -ne $n_index ]]
  then
     echo "ERROR: Bowtie2 indexing of genome $genome failed"
     exit 1
  fi
else
    echo "Skipping Bowtie2 indexing "
fi

cd  $tools_path/genomes/fasta
rm -f ${species_name}_mirna_fasta_u
rm -f mature.fa

