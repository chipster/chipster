#!/usr/bin/env bash

#
# This script will install Chipster 3, w/ dependencies
#
# Notice! This script needs super-user rights!!
# e.g. sudo bash install-chipster.sh 2>&1 | tee chipster.log
#

# Set execution trace
set -x

# Set exit on error
set -e

# Set fail on pipe
set -o pipefail

# Set exit on unset variable
set -u

# Function to do a hard retry with wget
wget_retry()
{
  while ! wget -t 1 -T 30 $*; do
    echo "wget failed, retrying"
    sleep 5
  done
}
export -f wget_retry

## System update
aptitude update
aptitude -y full-upgrade

## Install packages:

## Pre-requisites:

# Java
# sun-java6-jre (prefered would be openjdk-7-jre)
# Bypasses manual interaction to accept license
#echo 'sun-java6-jre shared/accepted-sun-dlj-v1-1 select true' | /usr/bin/debconf-set-selections
# aptitude -y install sun-java6-jre # (40 packages)
#aptitude -y --without-recommends install sun-java6-jre # (7 packages)
aptitude -y --without-recommends install openjdk-7-jre

# For now OpenJDK 6 is required, as is Ubuntu's default-jre,
# so we will set OpenJDK 7 as our selected choice
update-java-alternatives -s java-1.7.0-openjdk-amd64

# GNU Plot (w/o X)
# aptitude -y install gnuplot-nox # (55 packages)
aptitude -y --without-recommends install gnuplot-nox # (11 packages)

# Ghostscript
# Needed by (at least) qc-illumina.R
aptitude -y install ghostscript

## FASTX
aptitude -y install fastx-toolkit

## OpenMPI
#aptitude -y --without-recommends install openmpi-bin

## Python
# !! Anything from PyPI should/shall be installed with pip !!
# python-virtualenv
# virtualenvwrapper
aptitude -y --without-recommends install python-pip
# 2.6
# Needed for MACS?!?, IF REALLY THE CASE IT SHOULD BE RECOMPILED FOR 2.7
aptitude -y install python2.6

## Python Libraries:
# python-numpy, for HTSeq
# python-matplotlib, for HTSeq
aptitude -y --without-recommends install python-numpy python-matplotlib

## Libraries:
# build-essential (only devel)
# gfortran (libgfortran3)
# libcurl4-openssl-dev (libcurl3)
# libglib2.0-dev (libglib2.0-0)
# libglu1-mesa-dev (libglu1-mesa)
# libgsl0-dev (libgsl0ldbl)
# libpng-dev (libpng12-0)
# libreadline-dev (libreadline6) (>libreadline5-dev)
# libxml2-dev (libxml2)
# mesa-common-dev
# tcl-dev (tcl)
# tk-dev (tk)
# xorg-dev (only devel?)
# python-dev (python), for HTSeq
# libnetcdf6, for R
build_tools="yes" # Should tools be built, set to either "yes" or "no"
mode="devel" # Set to either "runtime" or "devel"

if [ $mode == "runtime" ]
then
  ## Runtime:
  aptitude -y --without-recommends install libgfortran3 libcurl3 libglib2.0-0 libglu1-mesa libgsl0ldbl libpng12-0 libreadline6 libxml2 mesa-common-dev tcl tk xorg-dev unixodbc gawk libnetcdf6 cython
elif [ $mode == "devel" ]
then
  ## Devel:
  aptitude -y --without-recommends install build-essential gfortran libcurl4-openssl-dev libglib2.0-dev libglu1-mesa-dev libgsl0-dev libpng-dev libreadline-dev libxml2-dev mesa-common-dev tcl-dev tk-dev xorg-dev python-dev unixodbc-dev libnetcdf-dev openjdk-7-jdk git ant cython
else
  echo "PROBLEM!!"
  exit 1
fi

## Perl Libraries:
# libjson-perl, for prinseq-graph
# libcairo-perl, for prinseq-graph
# libtext-simpletable-perl, for prinseq-graph
# libcontextual-return-perl, for prinseq-graph
# libwant-perl, for prinseq-graph
# cpanminus, for prinseq-graph
# Statistics::PCA, for prinseq-graph
# Math::Cephes, for prinseq-graph
# Math::MatrixReal, for prinseq-graph
aptitude -y --without-recommends install libjson-perl libcairo-perl libtext-simpletable-perl libcontextual-return-perl libwant-perl cpanminus
#cpanm Contextual::Return Exception::Class Test::{Warn,Exception,Differences,Deep} Math::Cephes Math::MatrixReal Statistics::PCA
cpanm --mirror http://ftp.funet.fi/pub/languages/perl/CPAN/ Contextual::Return Exception::Class Test::{Warn,Exception,Differences,Deep} Math::Cephes Math::MatrixReal Statistics::PCA

# misc packages
aptitude -y --without-recommends install unzip pigz pbzip2 dstat emacs23

## Initialize:
# Versions
CHIP_VER=3.0.0

# Paths
EXEC_PATH=${PWD}
INST_PATH=/opt
CHIP_PATH=${INST_PATH}/chipster
TOOLS_PATH=${CHIP_PATH}/tools
# check if TMPDIR is set and set it if necessary
if [ -z ${TMPDIR+x} ]; then
  TMPDIR=/tmp
fi
export TMPDIR
TMPDIR_PATH=$TMPDIR/install

# Misc
USERNAME=chipster
GROUPNAME=chipster

# nic.funet.fi service endpoint
NIC_MIRROR=bio.nic.funet.fi
#NIC_MIRROR=www.nic.funet.fi

# prebundled tar or git 
CHIPSTER_SOURCE="git"
#GIT_LABEL="tags/chipster-${CHIP_VER}" 
GIT_LABEL="data_beta" 

## Create tmpdir
rm -rf ${TMPDIR_PATH}
mkdir -p ${TMPDIR_PATH}/

## Chipster:

# Install Chipster

cd ${TMPDIR_PATH}/

# use prebuilt version
if [ $CHIPSTER_SOURCE == 'tar' ]; then
  
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/versions/${CHIP_VER}/chipster-${CHIP_VER}.tar.gz | tar -xz

elif [ $CHIPSTER_SOURCE == 'git' ]; then

  #build it from github
  rm -rf chipster
  #git clone https://github.com/chipster/chipster.git
  #cd chipster
  #git checkout $GIT_LABEL 
  curl -s -L http://github.com/chipster/chipster/tarball/$GIT_LABEL/ | tar -xz
  mv chipster-chipster-* chipster
  cd chipster
  
  rm -f keystore.ks
  keytool -genkey -alias chipster -keystore keystore.ks -storepass chipster -keypass chipster -validity 1825 -dname "cn=TEST, ou=TEST, o=TEST, c=TEST"
  echo "chipster" > passfeed; echo "chipster" >> passfeed; 
  ant < passfeed
  cd ..
  mv chipster/dist/chipster-${CHIP_VER}.tar.gz .
  rm -rf chipster
  tar xf chipster-${CHIP_VER}.tar.gz
fi

mv chipster/ ${CHIP_PATH}/

# Make some config "corrections"
sed -i'~' "s:/opt/chipster/:${CHIP_PATH}/:" ${CHIP_PATH}/comp/conf/runtimes.xml
# TODO The below should be made dynamic
sed -i'~' '/<configuration-module moduleId="comp">/a <!-- make compute service access filebroker file repository locally -->\
<entry entryKey="local-filebroker-user-data-path" type="string" \
description="path to local filebrokers user data directory">\
      <value>/opt/chipster/fileserver/file-root/user-data</value>\
</entry>' ${CHIP_PATH}/comp/conf/chipster-config.xml
sed -i'~' "s/#RUN_AS_USER=/RUN_AS_USER=${USERNAME}/" \
    ${CHIP_PATH}/activemq/bin/linux-x86-64/activemq \
    ${CHIP_PATH}/comp/bin/linux-x86-64/chipster-comp \
    ${CHIP_PATH}/auth/bin/linux-x86-64/chipster-auth \
    ${CHIP_PATH}/fileserver/bin/linux-x86-64/chipster-fileserver \
    ${CHIP_PATH}/webstart/bin/linux-x86-64/chipster-webstart \
    ${CHIP_PATH}/manager/bin/linux-x86-64/chipster-manager

# Make update.sh script available
cp ${CHIP_PATH}/admin/vm/update.sh ${CHIP_PATH}/update.sh
chmod u+x ${CHIP_PATH}/update.sh

# Symlink to tools
ln -s /mnt/tools ${TOOLS_PATH}

# Create user-data and jobs-data symlinks
mkdir ${CHIP_PATH}/fileserver/file-root/
ln -s /scratch/user-data ${CHIP_PATH}/fileserver/file-root/user-data
ln -s /scratch/jobs-data ${CHIP_PATH}/comp/jobs-data

# Symlink to genome browser annotations
mkdir ${CHIP_PATH}/fileserver/file-root/public/
ln -s ${TOOLS_PATH}/genomebrowser/annotations ${CHIP_PATH}/fileserver/file-root/public/annotations

touch ${CHIP_PATH}/auto-config-to-be-run

# TODO Chipster announcements
#mkdir /opt/chipster-announcements
#/opt/chipster-announcements/get-chipster-announcements.sh
# chown -R chipster:chipster /opt/chipster-announcements

#/etc/cron.d/chipster-announcements.crontab
#/etc/update-motd.d/31-vm-instructions
#/etc/update-motd.d/32-announcements
echo "Chipster announcements not done yet"

##############################################
# Install external applications and datasets #
##############################################


## In root:

# MACS, Artistic license
# part 1
cd ${TMPDIR_PATH}/
curl -s https://cloud.github.com/downloads/taoliu/MACS/MACS-1.4.2-1.tar.gz | tar -xz
cd MACS-1.4.2
python setup.py install
cd ..
rm -rf MACS-1.4.2
rm -f MACS-1.4.2-1.tar.gz

# MACS2, BSD licence
# part 1

cd $TMPDIR_PATH
curl -s 'https://cloud.github.com/downloads/taoliu/MACS/MACS-2.0.9-1.tar.gz' | tar -xz
cd MACS-2.0.9/
sudo python setup.py install
cd ..
rm -rf MACS-2.0.9

# HTSeq, GPL v3 or later
# part 1
pip install HTSeq==0.6.1
wget_retry -O /usr/local/bin/htseq-count_chr http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/htseq/htseq-count_chr 
chmod 755 /usr/local/bin/htseq-count_chr
wget_retry -O /usr/local/lib/python2.7/dist-packages/HTSeq/scripts/count_chr.py http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/htseq/count_chr_v2.py


## In tools:

if [ $mode == "devel" -a $build_tools == "yes" ]
then
  ## R
  
  # old R versions for Chipster 3.0

  # R-2.12.1-medips
  #curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.12.1-medips.tar.gz | tar -xz -C ${TOOLS_PATH}/

  # R-2.15.1-variantannotation
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1-variantannotation.tar.gz | tar -xz -C ${TOOLS_PATH}/

  
  R_VER=3.0.2  
  cd ${TMPDIR_PATH}/
  curl -s http://ftp.sunet.se/pub/lang/CRAN/src/base/R-3/R-${R_VER}.tar.gz | tar -xz
  cd R-${R_VER}/
  export MAKEFLAGS=-j
  ./configure --prefix=${TOOLS_PATH}/R-${R_VER}
  make
  make install
  echo 'MAKEFLAGS=-j' > ${TOOLS_PATH}/R-${R_VER}/lib/R/etc/Makevars.site # (could also be $HOME/.R/Makevars)
  # clean makeflags after R install
  export MAKEFLAGS=

  cd ../
  rm -rf R-${R_VER}/
  
  ${TOOLS_PATH}/R-${R_VER}/bin/Rscript --vanilla ${CHIP_PATH}/comp/modules/admin/R/install-libs.R

  # extra data for zinba R library
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/zinba-extras.tar.gz | tar xz -C ${TOOLS_PATH}

  ## R-2.15.1_bioc-2.11  
  #curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1_bioc-2.11/R-2.15.1_bioc-2.11-vmbin_v4.tar.gz | tar -xz -C ${TOOLS_PATH}/  

 	# Add RmiR.Hs.miRNA to R-2.15.1_bioc-2.11
  #curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1_bioc-2.11/library/RmiR.Hs.miRNA-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.15.1_bioc-2.11/lib64/R/library/

  ## External apps

  # Link tool admin scripts from Chipster installation
  mkdir ${TOOLS_PATH}/admin/
  ln -s ${CHIP_PATH}/comp/modules/ngs/admin ${TOOLS_PATH}/admin/ngs

  # Weeder, custom license, according to developers VM bundling is ok
  cd ${TMPDIR_PATH}/
  curl -s http://159.149.160.51/modtools/downloads/weeder1.4.2.tar.gz | tar -xz
  cd Weeder1.4.2/
  ./compileall
  cd ../
  mkdir ${TOOLS_PATH}/weeder/
  mv Weeder1.4.2/ ${TOOLS_PATH}/weeder/

  # ClusterBuster, no license
  cd ${TMPDIR_PATH}/
  wget_retry -O cbust-src.tar.gz http://zlab.bu.edu/~mfrith/downloads/cbust-src.tar.gz
  tar xf cbust-src.tar.gz
  rm -rf cbust-src.tar.gz
  cd cbust-src/
  make
  mkdir ${TOOLS_PATH}/ClusterBuster/
  mv cbust ${TOOLS_PATH}/ClusterBuster/
  cd ../
  rm -rf cbust-src/

  # Jaspar, no license
  cd ${TMPDIR_PATH}/
  wget_retry -nv http://zlab.bu.edu/clover/jaspar2
  mv jaspar2 ${TOOLS_PATH}/ClusterBuster/jaspar2005core.txt

  # Promoter sequence files, license?
  cd ${TMPDIR_PATH}/
  mkdir ${TOOLS_PATH}/weeder/seqs/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/weeder_seqs/All_Weeder_sequence_files_v1.tar.gz | tar -xz -C ${TOOLS_PATH}/weeder/seqs/

  # BEDTools, GNU GPL v2
  cd ${TMPDIR_PATH}/
  curl -s http://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz | tar -xz
  cd bedtools-2.17.0
  make clean
  make all
  cd ../
  mv bedtools-2.17.0 ${TOOLS_PATH}/
  ln -s bedtools-2.17.0 ${TOOLS_PATH}/bedtools

  #curl -L http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bedtools-2.17.0-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/


  # MACS 1 & 2 
  # part 2 - creating links
  cd ${TMPDIR_PATH}/
  mkdir ${TOOLS_PATH}/macs/
  ln -s /usr/local/bin/macs14 ${TOOLS_PATH}/macs/macs14
  ln -s /usr/local/bin/macs2 /mnt/tools/macs/macs2
  
  # SAM tools, BSD License, MIT License
  cd ${TMPDIR_PATH}/
  #curl -sL http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download | tar -xj
  #cd samtools-0.1.19/
  #make
  #cd ../
  #mv samtools-0.1.19/ ${TOOLS_PATH}
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/samtools-0.1.19-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s samtools-0.1.19 ${TOOLS_PATH}/samtools

  # tabix
  # curl -sL http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2/download | tar xj -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/tabix-0.2.6-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s tabix-0.2.6 ${TOOLS_PATH}/tabix 

  # Bowtie, Artistic License
  cd ${TMPDIR_PATH}/
  wget_retry -nv -O bowtie-0.12.7-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download
  unzip -q bowtie-0.12.7-linux-x86_64.zip
  mv bowtie-0.12.7/ ${TOOLS_PATH}
  ln -s bowtie-0.12.7 ${TOOLS_PATH}/bowtie
  rm bowtie-0.12.7-linux-x86_64.zip
  # remove example index
  rm ${TOOLS_PATH}/bowtie/indexes/e_coli.*

  # Bowtie 2, Artistic License
  cd ${TMPDIR_PATH}/
  wget_retry -nv -O bowtie2-2.1.0-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip/download
  unzip -q bowtie2-2.1.0-linux-x86_64.zip
  mv bowtie2-2.1.0 ${TOOLS_PATH}
  ln -s bowtie2-2.1.0 ${TOOLS_PATH}/bowtie2
  rm bowtie2-2.1.0-linux-x86_64.zip
  
  mkdir -p ${TOOLS_PATH}/bowtie2/indexes/
		
  # GRCh37_74 ensembl transcripts   
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/GRCh37_74_ensembl_transcripts.tar.gz | tar -xzv -C ${TOOLS_PATH}/bowtie2/indexes/
																																																														
  # FastQC, GPL v3 or later
  cd ${TMPDIR_PATH}/
  wget_retry -nv http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/fastqc_v0.10.0.zip
  unzip -q fastqc_v0.10.0.zip
  chmod a+x FastQC/fastqc
  mv FastQC/ ${TOOLS_PATH}/
  rm fastqc_v0.10.0.zip

  # HTSeq, GPL v3 or later
  # part 2
  cd ${TMPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/htseq/
  ln -s /usr/local/bin/htseq-qa ${TOOLS_PATH}/htseq/htseq-qa
  ln -s /usr/local/bin/htseq-count ${TOOLS_PATH}/htseq/htseq-count
  ln -s /usr/local/bin/htseq-count_chr ${TOOLS_PATH}/htseq/htseq-count_chr

  # Cufflinks, Boost License
  cd ${TMPDIR_PATH}/
  curl -s http://cufflinks.cbcb.umd.edu/downloads/cufflinks-1.0.3.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/  
  ln -s cufflinks-1.0.3.Linux_x86_64 ${TOOLS_PATH}/cufflinks
  curl -s http://cufflinks.cbcb.umd.edu/downloads/cufflinks-2.1.1.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s cufflinks-2.1.1.Linux_x86_64 ${TOOLS_PATH}/cufflinks2
 
  # Tophat, The Artistic License
  cd ${TMPDIR_PATH}/
  wget_retry -O tophat-1.3.2.Linux_x86_64.tar.gz http://ccb.jhu.edu/software/tophat/downloads/tophat-1.3.2.Linux_x86_64.tar.gz
  tar -xf tophat-1.3.2.Linux_x86_64.tar.gz
  rm -f tophat-1.3.2.Linux_x86_64.tar.gz
  mv tophat-1.3.2.Linux_x86_64 ${TOOLS_PATH}/
  ln -s tophat-1.3.2.Linux_x86_64 ${TOOLS_PATH}/tophat

  # Tophat 2, The Artistic License
  cd ${TMPDIR_PATH}/
  wget_retry -O tophat-2.0.10.Linux_x86_64.tar.gz http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.10.Linux_x86_64.tar.gz
  tar -xf tophat-2.0.10.Linux_x86_64.tar.gz -C ${TOOLS_PATH}/
  rm -f tophat-2.0.10.Linux_x86_64.tar.gz
  ln -s tophat-2.0.10.Linux_x86_64 ${TOOLS_PATH}/tophat2

  # BWA, GPL v3 or later, MIT License
  cd ${TMPDIR_PATH}/
  wget_retry -O bwa-0.6.1.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.6.1.tar.bz2/download
  tar -xf bwa-0.6.1.tar.bz2 --use-compress-program=pbzip2
  rm -f bwa-0.6.1.tar.bz2
  cd bwa-0.6.1/
  make
  cd ../
  mv bwa-0.6.1/ ${TOOLS_PATH}/
  ln -s bwa-0.6.1 ${TOOLS_PATH}/bwa
  
  # Fastx links
  mkdir -p ${TOOLS_PATH}/fastx/bin/
  ln -s /usr/bin/fasta_* ${TOOLS_PATH}/fastx/bin/
  ln -s /usr/bin/fastq_* ${TOOLS_PATH}/fastx/bin/
  ln -s /usr/bin/fastx_* ${TOOLS_PATH}/fastx/bin/

  # prinseq
  cd ${TMPDIR_PATH}/
  curl -sL http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz/download | tar -xz
  chmod a+x prinseq-lite-0.20.4/*.pl
  mv prinseq-lite-0.20.4 ${TOOLS_PATH}/
  ln -s prinseq-lite-0.20.4 ${TOOLS_PATH}/prinseq

  # Genome data for tools
  cd ${TMPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/genomes/

  # GTF gene data for tools
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Rattus_norvegicus.RGSC3.4.68.tar.gz | tar -xz -C ${TOOLS_PATH}/

  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Canis_familiaris.CanFam3.1.71.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Gasterosteus_aculeatus.BROADS1.71.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar xz -C ${TOOLS_PATH}/

  # miRNA mapping data
  cd ${TMPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/miRNA_mappings/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/miRNA_mappings/All_miRNA_mappings_v1.tar.gz | tar -xz -C ${TOOLS_PATH}/miRNA_mappings/

  # Bowtie indexes, built for Chipster
  cd ${TMPDIR_PATH}/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_athaliana.TAIR10.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_mm9.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_miRBase19_Rattus_norvegicus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_miRBase19_Homo_sapiens.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_miRBase19_Mus_musculus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_Canis_familiaris.CanFam3.1.71.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_Gasterosteus_aculeatus.BROADS1.71.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/

  ln -s -t ${TOOLS_PATH}/bowtie/indexes ../../genomes/fasta/nochr
  rm -f ${TOOLS_PATH}/bowtie/indexes/Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa
  ln -s ../../genomes/fasta/nochr/Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa ${TOOLS_PATH}/bowtie/indexes

	# Bowtie2 indexes, built for Chipster
  cd ${TMPDIR_PATH}/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_athaliana.TAIR10.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_mm9.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_miRBase19_Rattus_norvegicus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_miRBase19_Homo_sapiens.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_miRBase19_Mus_musculus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/	
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Canis_familiaris.CanFam3.1.71.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Gasterosteus_aculeatus.BROADS1.71.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
	
	ln -s -t ${TOOLS_PATH}/bowtie2/indexes ../../genomes/fasta/nochr
  rm -f ${TOOLS_PATH}/bowtie2/indexes/Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa
  ln -s ../../genomes/fasta/nochr/Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa ${TOOLS_PATH}/bowtie2/indexes

  # Tophat indexes
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/hg19.ti.tar.gz | tar -xzv -C ${TOOLS_PATH}/bowtie2/indexes/

  # bwa indexes, built for Chipster
  cd ${TMPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/bwa_indexes/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_mm9.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_athaliana.TAIR10.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_Canis_familiaris.CanFam3.1.71.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_Gasterosteus_aculeatus.BROADS1.71.tar.gz | tar xz -C ${TOOLS_PATH}/

  # Fasta files for tools

  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_mm9.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_athaliana.TAIR10.tar.gz | tar xz -C ${TOOLS_PATH}/

  # Fasta files for genome browser

  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Arabidopsis_thaliana.TAIR10.17.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Bos_taurus.UMD3.1.70.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Canis_familiaris.CanFam3.1.70.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Canis_familiaris.CanFam3.1.71.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Gallus_gallus.Gallus_gallus-4.0.pre-v2.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Gasterosteus_aculeatus.BROADS1.70.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Gasterosteus_aculeatus.BROADS1.71.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Vitis_vinifera.IGGP_12x.17.tar.gz | tar -xz -C ${TOOLS_PATH}/

  # Genome annotations for genome browser
  cd ${TMPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/genomebrowser/annotations/

  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/annotations/compressed/3.0.0/Arabidopsis_thaliana.TAIR10.17.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/annotations/compressed/3.0.0/Bos_taurus.UMD3.1.70.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/annotations/compressed/3.0.0/Canis_familiaris.CanFam3.1.70.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/annotations/compressed/3.0.0/Gallus_gallus.Gallus_gallus-4.0.pre.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/annotations/compressed/3.0.0/Gasterosteus_aculeatus.BROADS1.70.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/annotations/compressed/3.0.0/Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/annotations/compressed/3.0.0/Vitis_vinifera.IGGP_12x.17.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/

    # Genome bundles
  cd ${CHIP_PATH}/
  apt-get -y install python3-yaml #sudo
  wget_retry http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/bundle/bundles.yaml -O bundles.yaml
  #python3 bundle.py install all
  python3 bundle.py install Drosophila_melanogaster.BDGP5.bowtie
  python3 bundle.py install Drosophila_melanogaster.BDGP5.bowtie2
  python3 bundle.py install Drosophila_melanogaster.BDGP5.bwa
  python3 bundle.py install Drosophila_melanogaster.BDGP5.gb
  python3 bundle.py install Drosophila_melanogaster.BDGP5

  python3 bundle.py install Schizosaccharomyces_pombe.ASM294v2.bowtie
  python3 bundle.py install Schizosaccharomyces_pombe.ASM294v2.bowtie2
  python3 bundle.py install Schizosaccharomyces_pombe.ASM294v2.bwa
  python3 bundle.py install Schizosaccharomyces_pombe.ASM294v2.gb
  python3 bundle.py install Schizosaccharomyces_pombe.ASM294v2

  python3 bundle.py install Canis_familiaris.BROADD2
  python3 bundle.py install Escherichia_coli_n1.GCA_000303635.1.18
  python3 bundle.py install Halorubrum_lacusprofundi_ATCC_49239
  python3 bundle.py install Homo_sapiens.GRCh37
  python3 bundle.py install Homo_sapiens.NCBI36
  python3 bundle.py install Human-MT.NC_012920.1
  python3 bundle.py install Mus_musculus.GRCm38
  python3 bundle.py install Mus_musculus.NCBIM37
  python3 bundle.py install Ovis_aries.Oar_v3.1
  python3 bundle.py install Rattus_norvegicus.RGSC3.4
  python3 bundle.py install Sus_scrofa.Sscrofa10.2

  # Some bundles are missing few files

  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Human-MT.NC_012920.1-v2.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Mus_musculus.GRCm38.68.tar.gz | tar xz -C ${TOOLS_PATH}/ genomes/gtf/Mus_musculus.NCBIM37.62.chr.gtf
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Homo_sapiens.GRCh37.68.tar.gz | tar xz -C ${TOOLS_PATH}/ genomes/gtf/Homo_sapiens.GRCh37.68.DEXSeq.gtf genomes/gtf/Homo_sapiens.GRCh37.68.chr.DEXSeq.gtf


  # ConsensuPathDB
  apt-get -y install python-zsi #sudo

  # DEXSeq
	cd ${TMPDIR_PATH}/
	#	curl -sL http://www.bioconductor.org/packages/release/bioc/src/contrib/DEXSeq_1.8.0.tar.gz | tar -xz
	curl -sL http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/DEXSeq_1.8.0.tar.gz | tar -xz
	mkdir ${TOOLS_PATH}/dexseq-exoncounts
	cp DEXSeq/inst/python_scripts/dexseq_count.py ${TOOLS_PATH}/dexseq-exoncounts  
	cp DEXSeq/inst/python_scripts/dexseq_prepare_annotation.py ${TOOLS_PATH}/dexseq-exoncounts
	rm -rf DEXSeq	

  # vcftools, GPLv3
  cd ${TMPDIR_PATH}/
  wget_retry -O vcftools_0.1.11.tar.gz http://sourceforge.net/projects/vcftools/files/vcftools_0.1.11.tar.gz/download
  tar xf vcftools_0.1.11.tar.gz
  cd vcftools_0.1.11/
  make
  cd ../
  rm -f vcftools_0.1.11.tar.gz
  mv vcftools_0.1.11/ ${TOOLS_PATH}/
  #curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/vcftools_0.1.11-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s vcftools_0.1.11 ${TOOLS_PATH}/vcftools


  # GATK, MIT
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/GenomeAnalysisTKLite-latest.tar.bz2 | tar -xj -C ${TOOLS_PATH}/
  ln -s GenomeAnalysisTKLite-2.1-11-gfb37f33 ${TOOLS_PATH}/GenomeAnalysisTK2

  # tagcleaner, GPLv3
  cd ${TMPDIR_PATH}/
  wget_retry -nv -O tagcleaner-standalone-0.12.tar.gz http://downloads.sourceforge.net/project/tagcleaner/standalone/tagcleaner-standalone-0.12.tar.gz
  tar xf tagcleaner-standalone-0.12.tar.gz -C ${TOOLS_PATH}/
  rm -f tagcleaner-standalone-0.12.tar.gz
  ln -s tagcleaner-standalone-0.12 ${TOOLS_PATH}/tagcleaner

  # fseq, GPLv3
  curl -s http://fureylab.med.unc.edu/fseq/fseq_1.84.tgz | tar -xz -C ${TMPDIR_PATH}/ 
  mv ${TMPDIR_PATH}/fseq ${TOOLS_PATH}/fseq-1.84
  ln -s fseq-1.84 ${TOOLS_PATH}/fseq
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/read_extend_bed.pm -o ${TOOLS_PATH}/fseq/bin/read_extend_bed.pm
  chmod 775 ${TOOLS_PATH}/fseq/bin/read_extend_bed.pm

  # mothur GPLv3
  cd ${TMPDIR_PATH}/
  wget_retry -nv http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/mothur/Mothur-1.28.cen_64.noReadLine.zip
  unzip -q Mothur-1.28.cen_64.noReadLine.zip
  mv mothur ${TOOLS_PATH}/mothur-1.28
  ln -s mothur-1.28 ${TOOLS_PATH}/mothur
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/mothur/mothur-data.tar.gz | tar -xz -C ${TOOLS_PATH}/

  # Picard tools, Apache License V2.0, MIT
  cd ${TMPDIR_PATH}/
  wget_retry -nv -O picard-tools-1.105.zip http://sourceforge.net/projects/picard/files/picard-tools/1.105/picard-tools-1.105.zip/download
  unzip -q picard-tools-1.105.zip
  rm picard-tools-1.105.zip
  # remove this optional jar because it's in the root of the zip
  rm snappy-java-1.0.3-rc3.jar
  mv picard-tools-1.105/ ${TOOLS_PATH}
  cd ${TOOLS_PATH}
  ln -s picard-tools-1.105 picard-tools

  # RSeQC, GPLv3
  cd ${TMPDIR_PATH}/
  wget_retry -nv -O RSeQC-2.3.7.tar.gz http://sourceforge.net/projects/rseqc/files/RSeQC-2.3.7.tar.gz/download
  tar xf RSeQC-2.3.7.tar.gz
  mv RSeQC-2.3.7 ${TOOLS_PATH}/RSeQC-2.3.7
  rm -f RSeQC-2.3.7.tar.gz
  cd ${TOOLS_PATH}
  rm -f RSeQC
  ln -s RSeQC-2.3.7 RSeQC
  cd RSeQC
  python setup.py install #sudo

   # EMBOSS, GPL
  apt-get -y install libgd2-noxpm-dev # sudo, emboss needs this to create png images
  # also vmbin from nic

#  EMBOSS_MIRROR=ftp://emboss.open-bio.org/pub/EMBOSS/old/6.5.0
  EMBOSS_MIRROR=http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS

  EMBOSS_VERSION=6.5.7
  EMBOSS_PATH=${TOOLS_PATH}/EMBOSS-${EMBOSS_VERSION}
  # note version in path                                                                                                                                                              
	curl ${EMBOSS_MIRROR}/EMBOSS-${EMBOSS_VERSION}.tar.gz | tar -xz -C ${TMPDIR_PATH}/
  cd ${TMPDIR_PATH}/EMBOSS-6.5.7
	
  #wget ftp://emboss.open-bio.org/pub/EMBOSS/fixes/patches/patch-1-11.gz                                                                                                              
  #gunzip patch-1-11.gz                                                                                                                                                               
  #patch -p1 < patch-1-11                                                                                                                                                             
	
	EMBOSS_OPTIONS="--prefix=${EMBOSS_PATH}"
	./configure ${EMBOSS_OPTIONS}
	make
	make install
  ln -s EMBOSS-6.5.7 ${TOOLS_PATH}/emboss


  # EMBOSS extras

	curl ${EMBOSS_MIRROR}/MEME-4.7.650.tar.gz | tar -xz -C ${TMPDIR_PATH}/
	cd ${TMPDIR_PATH}/MEME-4.7.650
	./configure ${EMBOSS_OPTIONS}
	make
	make install
	cd ..

  # phylipnew                                                                                                                                                                         
	curl ${EMBOSS_MIRROR}/PHYLIPNEW-3.69.650.tar.gz | tar -xz -C ${TMPDIR_PATH}/
	cd PHYLIPNEW-3.69.650
	./configure ${EMBOSS_OPTIONS}
	make
	make install
	cd ..
	
  # vienna                                                                                                                                                                            
	curl ${EMBOSS_MIRROR}/VIENNA-1.7.2.650.tar.gz | tar -xz -C ${TMPDIR_PATH}/
	cd  ${TMPDIR_PATH}/VIENNA-1.7.2.650
	./configure ${EMBOSS_OPTIONS}
	make
	make install
	cd ..

	# REBASE reference data and indeces	
	cd ${EMBOSS_PATH}/share/EMBOSS/data/REBASE
	wget_retry ftp://ftp.neb.com/pub/rebase/withrefm.txt
	wget_retry ftp://ftp.neb.com/pub/rebase/proto.txt
	../../../../bin/rebaseextract -infile withrefm.txt -protofile proto.txt
	
  # primer3, BSD                                                                                                                                                                      
	mkdir ${TOOLS_PATH}/primer3
	curl -L http://sourceforge.net/projects/primer3/files/primer3/1.1.4/primer3-1.1.4.tar.gz | tar -xz -C ${TOOLS_PATH}/primer3
	cd ${TOOLS_PATH}/primer3/src/
	make
	
  # meme, quite free, see documentation                                                                                                                                               
	cd ${TMPDIR_PATH}
	curl http://ebi.edu.au/ftp/software/MEME/4.2.0/meme_4.2.0.tar.gz | tar -xz -C ${TMPDIR_PATH}
	cd meme_4.2.0/
	./configure --prefix=${TOOLS_PATH}/meme_4.2.0 --with-url="http://meme.nbcr.net/meme"
	make
	make install
	ln -s meme_4.2.0 ${TOOLS_PATH}/meme
	cd ..
	rm -rf ${TMPDIR_PATH}/meme_4.2.0

  # dimont, GPLv3?
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/dimont.tar.gz | tar -xz -C ${TOOLS_PATH}/

  # Trimmomatic, GPL
  cd ${TMPDIR_PATH}/
  wget_retry -nv http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
  unzip Trimmomatic-0.32.zip
  mv Trimmomatic-0.32 ${TOOLS_PATH}/
  ln -s Trimmomatic-0.32 ${TOOLS_PATH}/trimmomatic
  rm Trimmomatic-0.32.zip
  
  # Express, Artistic license 2.0
  curl http://bio.math.berkeley.edu/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz | tar -xz -C ${TOOLS_PATH}/
  ln -s express-1.5.1-linux_x86_64 ${TOOLS_PATH}/express
  
  # BLAST, public domain
  curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.29+-x64-linux.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s ncbi-blast-2.2.29+ ${TOOLS_PATH}/blast                                  

  # Mafft, BSD
  #cd ${TMPDIR_PATH}/
  #curl http://mafft.cbrc.jp/alignment/software/mafft-7.130-without-extensions-src.tgz | tar -xz
  #cd mafft-7.130-without-extensions/core
  #sed -i 's/PREFIX = \/usr\/local/PREFIX = \/opt\/chipster\/tools\/mafft-7.130-without-extensions/g' Makefile
  #make clean
  #make
  #make install
  curl http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/mafft-7.130-without-extensions-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s mafft-7.130-without-extensions ${TOOLS_PATH}/mafft

  # QDNAseq files from Ilari
  curl -s http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/misc/QDNAseq.hg19.tar.gz | tar -xz -C ${TOOLS_PATH}/



  ## Create checksums
  cd ${TOOLS_PATH}/
  find . '!' -type d '!' -type l -print0 | xargs -0 sha256sum >> tools.sha256sum

fi

## Fix rights 
chown -R -h ${USERNAME}.${GROUPNAME} ${CHIP_PATH}/
if [ -d /mnt/tools/ ]; then
  chown -R ${USERNAME}.${GROUPNAME} /mnt/tools/
  chmod -R go-rwxst,u=rwX,go=rX /mnt/tools/
fi

## Clean up:
ls -lah ${TMPDIR_PATH}/ # Debug, list uncleaned mess
rm -rf ${TMPDIR_PATH}/

## Init.d:
ln -s ${CHIP_PATH}/chipster /etc/init.d/chipster
update-rc.d chipster defaults 30 70 # start 99 3 5 . stop 1 0 1 2 4 6 .

ln -s ${CHIP_PATH}/activemq/bin/linux-x86-64/activemq /etc/init.d/chipster-activemq
ln -s ${CHIP_PATH}/comp/bin/linux-x86-64/chipster-comp /etc/init.d/chipster-comp
ln -s ${CHIP_PATH}/auth/bin/linux-x86-64/chipster-auth /etc/init.d/chipster-auth
ln -s ${CHIP_PATH}/fileserver/bin/linux-x86-64/chipster-fileserver /etc/init.d/chipster-fileserver
ln -s ${CHIP_PATH}/webstart/bin/linux-x86-64/chipster-webstart /etc/init.d/chipster-webstart
ln -s ${CHIP_PATH}/manager/bin/linux-x86-64/chipster-manager /etc/init.d/chipster-manager

