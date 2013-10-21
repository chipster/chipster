#!/bin/bash
# This script packages genomes to bundle packages. First, use script 'add_genome.sh' to fetch and 
# convert necessary files. Then check that all files are ok and fill in the details of each genome 
# in the yaml file under tools/genomebrowser/annotation. After that, this script can be used to create the bundle
# packages. Internally the packages are created with a tool 'to_bundle.py'.

cd /opt/chipster/tools

cat genomes/genome_list | cut -d " " -f 2,4 | sort > genomes.txt

#Some old genomes
#echo "Canis_familiaris BROADD2.67" >> genomes.txt
#echo "Homo_sapiens NCBI36.54" >> genomes.txt
#echo "Mus_musculus NCBIM37.67" >> genomes.txt
#echo "Rattus_norvegicus RGSC3.4.69" >> genomes.txt

#echo "Human_mitoch NC_012920" >> genomes.txt
#echo "Ovis_aries Oar_v3.1" >> genomes.txt
#echo "phiR1-RT HE956709.1" >> genomes.txt
#echo "Yersinia_similis N916Ysi" >> genomes.txt

BUNDLE_VERSION="1.1"
CHIPSTER_VERSION="2.9"

mkdir bundle-file-lists

find /opt/chipster/tools/* -type f > all.txt

while read LINE; do
  SPECIES=$(echo $LINE | cut -d " " -f 1)
  VERSION=$(echo $LINE | cut -d " " -f 2)

  #file list names
  BOWTIE="bundle-file-lists/$SPECIES.$VERSION.bowtie.txt"
  BOWTIE2="bundle-file-lists/$SPECIES.$VERSION.bowtie2.txt"
  BWA="bundle-file-lists/$SPECIES.$VERSION.bwa.txt"
  OTHER="bundle-file-lists/$SPECIES.$VERSION.txt"

  #create file lists
  #cat all.txt | grep $SPECIES | grep $VERSION | grep "bowtie-" > $BOWTIE
  #cat all.txt | grep $SPECIES | grep $VERSION | grep "bowtie2" > $BOWTIE2
  #cat all.txt | grep $SPECIES | grep $VERSION | grep "bwa_indexes" > $BWA
  #cat all.txt | grep $SPECIES | grep $VERSION | grep -v "bowtie-" | grep -v "bowtie2" | grep -v "bwa_indexes" >  $OTHER

  #create bundles
  #if file list is not empty
  if [ -s $BOWTIE ]
  then
    : #empty script block when following is commented out
    cat $BOWTIE | python3 to_bundle.py -n $SPECIES.$VERSION.bowtie -v $BUNDLE_VERSION -p $CHIPSTER_VERSION
  fi

  if [ -s $BOWTIE2 ]
  then
    : #empty script block when following is commented out
    cat $BOWTIE2 | python3 to_bundle.py -n $SPECIES.$VERSION.bowtie2 -v $BUNDLE_VERSION -p $CHIPSTER_VERSION
  fi

  if [ -s $BWA ]
  then
    : #empty script block when following is commented out
    cat $BWA | python3 to_bundle.py -n $SPECIES.$VERSION.bwa -v $BUNDLE_VERSION -p $CHIPSTER_VERSION
  fi

  if [ -s $OTHER ]
  then
    : #empty script block when following is commented out
    cat $OTHER | python3 to_bundle.py -n $SPECIES.$VERSION -v $BUNDLE_VERSION -p $CHIPSTER_VERSION
  fi

done < genomes.txt


#Some legacy genomes have non-standard names. To create bundles for these, 
#you have to run following command  after the file lists are created, but before to_bundle is run

#find /opt/chipster/tools/* | grep mm9 | grep fasta >> bundle-file-lists/Mus_musculus.NCBIM37.67.txt
#find /opt/chipster/tools/* | grep rn4 | grep fasta >> bundle-file-lists/Rattus_norvegicus.RGSC3.4.69.txt 
#find /opt/chipster/tools/* | grep Human-MT | grep fasta >> bundle-file-lists/Human_mitoch.NC_012920.txt
#find /opt/chipster/tools/* | grep ovis_aries | grep fasta >> bundle-file-lists/Ovis_aries.Oar_v3.1.txt 
#find /opt/chipster/tools/* | grep R1-RT | grep fasta >> bundle-file-lists/phiR1-RT.HE956709.1.txt
#find /opt/chipster/tools/* | grep N916Ysi | grep fasta >> bundle-file-lists/Yersinia_similis.N916Ysi.txt 

rm genomes.txt
rm all.txt

# to_bundle.py creates separate yaml files
cat *.yaml > bundles-local.yaml
# prefix package file names with url, sed explanation:
# - search for 'packages:'
# - 'n' to continue processing on next line
# - 's' to search and replace
# - search for first four spaces 
# - replace with four spaces and url
sed '/packages:/{n; s/    /    http:\/\/www.nic.funet.fi\/pub\/sci\/molbio\/chipster\/dist\/tools_extras\/bundle\//}' bundles-local.yaml > bundles-$BUNDLE_VERSION.yaml

#If everything went well
# - copy packages and bundles-X.X.yaml to nic /pub/sci/molbio/chipster/dist/tools_extras/bundle/
# - append bundles-X.X.yaml to bundles.yaml


