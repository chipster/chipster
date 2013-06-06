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


while [[ $# -ge 1 ]]
do
  case "$1" in
              '-species')
                species="$2"
                ensembl=0
                shift 
                shift 
              ;;
              '-fasta')
                genome_fasta="$2"
                fasta=1
                shift 
                shift 
              ;;
              '-gtf')
                genome_gtf="$2"
                gtf=1
                shift 
                shift 
              ;;
              '-version')
                genome_version="$2"
                shift 
                shift 
              ;;
              '-descfilter')
                desctext="$2"
                text="1"
                shift 
                shift 
              ;;
              '-minlength')
                length=1
                minlength="$2"
                shift
                shift
              ;;
              '-only_bwa')
                INDEX_BOWTIE=0
                INDEX_BOWTIE2=0
                shift
              ;;
	      '-only_bowtie2')
                INDEX_BOWTIE=0
                INDEX_BWA=0
                shift
              ;;

              *)
                species="$1"
                ensembl=1
                shift 
              ;; 

    esac
done


##
#Retrieve the fasta file
##

#Check if  taxid is used in stead of name

taxid=$(echo $species | tr -d "[a-z,A-Z]" )
species=$(echo $species | sed s/" "/"_"/g )
#test for taxnumber
if [[ "$species" == "$taxid" ]]
then
  species=$(taxget taxon:$taxid -oformat excel -filter | sed s/" "/"_"/g | awk '{print $5}')
  echo "Taxid: $taxid corresponds species: $species"
else
  tax_name=$(echo $species | sed s/"_"/" "/g )
  taxid=$(grep -i "|.$tax_name.|" $tools_path/emboss/share/EMBOSS/data/TAXONOMY/names.dmp  | awk '{print $1}')
fi


echo $species $taxid

#reading the data from ensembl

echo $ensembl

if [[ $ensembl -eq 1 ]]
then
  echo "Retrtieving and indexing genome sequence for $species"

  cd ${tools_path}/genomes/fasta/nochr
  echo ensemblfetch.sh $species
  genome_fasta=$(ensemblfetch.sh $species | tail -1)

  if [[ $genome_fasta == "--------------------------------------------------------------------------------" ]]
  then
    echo Species $species was not found from the Ensembl database.
    exit 1
  fi

  genome_release=$(echo $genome_fasta | awk -F "." '{print $2}')

  #Try to find the build name
  wget "http://hgdownload.soe.ucsc.edu/admin/hgcentral.sql"
  # Take only inserts of table dbDB
  cat hgcentral.sql |grep "INSERT INTO \`dbDb\`" > dbDb.txt
  # Take only specific columns
  cat dbDb.txt | cut -d "(" -f 4- |cut -d "'" -f 1,11,16 > release-ucsc.txt
  # Remove some characters
  sed -i 's/)//g' release-ucsc.txt
  sed -i 's/;//g' release-ucsc.txt
  sed -i 's/,//g' release-ucsc.txt
  # Replace slash and quotes with tab
  sed -i 's/\//	/g' release-ucsc.txt #note tab
  sed -i "s/'/	/g" release-ucsc.txt #note tab


  version=$(awk '{if ( $1=="'$genome_release'" ) if ( $NF=="'$taxid'") print $2}' release-ucsc.txt) 
  rm -f release-ucsc.txt
  rm -f hgcentral.sql
  rm -f dbDb.txt
  if [[ $version == "" ]]
  then
    version=$genome_release
  fi

fi



if [[ $fasta -eq 1 ]]
then
  cp $genome_fasta ${tools_path}/genomes/fasta/
  cd ${tools_path}/genomes/fasta
fi


##
#Check if the fasta file as already been indexed
##

size=$(ls -l $genome_fasta | awk '{print $5}')
checksum=$(md5sum $genome_fasta | awk '{print $1}')

#look for matching size and md5sum

genome_check=$(grep -h $size $tools_path/genomes/genome_list | grep $checksum | awk '{print $1}' | tail -1)


#
if [ ! $genome_check == "" ]; then
  echo "File $genome_fasta has alredy been indexed"
  exit 0
fi


if [[ "$text" == "1" ]]
  then
  textsearch  -sequence $genome_fasta -pattern "$desctext" -outfile outfile_"$$"_list -only -usa -auto >> /dev/null
  test=$(head outfile_"$$"_list | wc -l) 
  if [[ $test -eq 0 ]] 
  then 
     echo "No sequences matching search pattern: $desctext were found" > $outfile
     rm -rf outfile_"$$"_list
     exit 1
  fi
  seqret @outfile_"$$"_list -outseq ${genome_fasta}.filtered -auto
  rm -f $genome_fasta
  mv ${genome_fasta}.filtered $genome_fasta
fi



###
#  get the gtf file
###
if [[ $ensembl -eq 1 ]]
then
  cd ${tools_path}/genomes/gtf
  ensemblfetch.sh -type gtf $species
fi

if [[ $gtf -eq 1 ]]
then
  cp $location/$genome_gtf ${tools_path}/genomes/gtf/
fi



###

#make bwa_indexes
if [[ $INDEX_BWA -eq 1 ]]
then
  echo -------------------------------------------------------
  echo Calculating BWA indexes for $genome_fasta
  cd $tools_path/bwa_indexes
  if [[ $ensembl -eq 1 ]]
  then
    ln -s $tools_path/genomes/fasta/nochr/$genome_fasta $genome_fasta
  else 
    ln -s $tools_path/genomes/fasta/$genome_fasta $genome_fasta
  fi
  $tools_path/bwa/bwa index $genome_fasta
else
  echo "Skipping BWA indexing"
fi


#make bowtie_indexes
if [[ $INDEX_BOWTIE -eq 1 ]]
then
  echo -------------------------------------------------------
  echo Calculating Bowtie indexes for $genome_fasta 
  cd $tools_path/bowtie/indexes
  if [[ $ensembl -eq 1 ]]
  then
   ln -s $tools_path/genomes/fasta/nochr/$genome_fasta $genome_fasta
  else
   ln -s $tools_path/genomes/fasta/$genome_fasta $genome_fasta
  fi
  bowtie_name=$(basename $genome_fasta .fa)
  $tools_path/bowtie/bowtie-build $genome_fasta $bowtie_name

  #check bowtie2 indexes
  n_genome=$(grep -c "^>" $genome_fasta)
  n_index=$($tools_path/bowtie/bowtie-inspect -n $bowtie_name | wc -l)

  if [[ $n_genome -ne $n_index ]]
  then
     echo "ERROR: Bowtie indexing of genome $bowtie2_name failed"
     exit 1
  else
    echo Bowtie index of genome $bowtie_name OK
  fi
else
    echo "Skipping Bowtie indexing"
fi





#make bowtie2_indexes
if [[ $INDEX_BOWTIE2 -eq 1 ]]
then
  echo -------------------------------------------------------
  echo Calculating Bowtie2 indexes for $genome_fasta 

  cd $tools_path/bowtie2/indexes
  if [[ $ensembl -eq 1 ]]
  then
   ln -s $tools_path/genomes/fasta/nochr/$genome_fasta $genome_fasta
  else
    ln -s $tools_path/genomes/fasta/$genome_fasta $genome_fasta
  fi
  bowtie2_name=$(basename $genome_fasta .fa)
  $tools_path/bowtie2/bowtie2-build $genome_fasta $bowtie2_name

  #check bowtie2 indexes
  n_genome=$(grep -c "^>" $genome_fasta)
  n_index=$($tools_path/bowtie2/bowtie2-inspect -n $bowtie2_name | wc -l)

  if [[ $n_genome -ne $n_index ]]
  then
     echo "ERROR: Bowtie2 indexing of genome $bowtie2_name failed"
     exit 1
  else
    echo Bowtie2 index of genome $bowtie2_name OK
  fi
else
    echo "Skipping Bowtie indexing"
fi

day=$(date)
echo $taxid $species $genome_fasta $version $size $day $checksum >> $tools_path/genomes/genome_list
