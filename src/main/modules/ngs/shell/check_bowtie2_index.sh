#!/bin/bash 
#
#Automatic index cheking and indexing toold for bwa indexes
#NOTE! Does not do color space indexes
#K.M. 12.1. 2012

bowtie2_path=("/opt/chipster/tools/bowtie2")


if [ -d "/scratch/tmp" ]; then
   mkdir -p "/scratch/tmp/bowtie2_indexes/tmp"
   index_path=("/scratch/tmp/bowtie2_indexes/tmp")
else
   index_path=("./")
fi

genome=($1)
size=(`ls -l $genome | awk '{print $5}' `)
checksum=(`md5sum $genome | awk '{print $1}'`)
location=(`pwd`)

if [ ! -d $index_path/genome_0 ]; then
  mkdir $index_path/genome_0
fi


#look for matching size and md5sum
genome_dir=(`grep -h $size $index_path/genome_*/size_and_md5 | grep $checksum | awk '{print $1}' | tail -1`)


#
if [ ! $genome_dir == "" ]; then
   echo "Pre-indexed genome found"
else
   echo "Calculating indexes"
   cd $index_path
   genome_dir=(`ls -d genome_* | awk -F "_" '{ print $NF }' | sort -n | tail -1 | awk '{ a = ( $1 + 1 ) }{print "genome_"a}'`)
   mkdir $genome_dir
   cd $genome_dir
   cp "$location"/"$genome" ./$genome
  

   $bowtie2_path/bowtie2-build $genome $genome 


   #check that idexes are found
   echo "$genome"
   for f in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
       if [ -e "$genome"."$f" ]; then
         echo "$genome.$f OK"
       else
         echo "Indexing failed"
         echo "Index file $genome.$f not found"
         exit 1
       fi
   done
   echo "$genome_dir $size $checksum" > size_and_md5
   cd $location
fi

echo "The bowtie2_indexes are in directory:"
echo "$index_path/$genome_dir"
exit 

