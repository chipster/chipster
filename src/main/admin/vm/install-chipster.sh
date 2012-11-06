#!/usr/bin/env bash

#
# This script will install Chipster 2, w/ dependencies
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

# Wrapper function to run commands and check return values
#run ()
#{
#    # NOTICE!! Nothing is allowed to be sent to standard out!!
#
#    # Uncomment to enable debugging
#    debug=true
#    
#    # Print command
#    [ debug ] && echo "Run: Cmdline: $@" > &2
#    
#    # Run command
#    eval $@
#    retval=$?
#    
#    # Print return value
#    [ debug ] && echo "Run: Return value: ${retval}" > &2
#    
#    # Check return value
#    if [ ${retval} -ne 0 ]
#    then
#        [ debug ] && echo "Run: Command failed!"  > &2
#    else
#        [ debug ] && echo "Run: Command succeeded!"  > &2
#    fi
#}

## System update
aptitude update
aptitude -y full-upgrade

## Install packages:

## Pre-requesites:

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
build_tools="yes" # Should tools be built, set to either "yes" or "no"
mode="devel" # Set to either "runtime" or "devel"
if [ $mode == "runtime" ]
then
  ## Runtime:
  # aptitude -y install libgfortran3 libcurl3 libglib2.0-0 libglu1-mesa libgsl0ldbl libpng12-0 libreadline6 libxml2 mesa-common-dev tcl tk xorg-dev (141 packages)
  aptitude -y --without-recommends install libgfortran3 libcurl3 libglib2.0-0 libglu1-mesa libgsl0ldbl libpng12-0 libreadline6 libxml2 mesa-common-dev tcl tk xorg-dev unixodbc gawk # (117 packages)
elif [ $mode == "devel" ]
then
  ## Devel:
  # aptitude -y install build-essential gfortran libcurl4-openssl-dev libglib2.0-dev libglu1-mesa-dev libgsl0-dev libpng-dev libreadline-dev libxml2-dev mesa-common-dev tcl-dev tk-dev xorg-dev (181 packages)
  aptitude -y --without-recommends install build-essential gfortran libcurl4-openssl-dev libglib2.0-dev libglu1-mesa-dev libgsl0-dev libpng-dev libreadline-dev libxml2-dev mesa-common-dev tcl-dev tk-dev xorg-dev python-dev unixodbc-dev #libopenmpi-dev # (167 packages)
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
cpanm Statistics::PCA Math::Cephes Math::MatrixReal 

## Initialize:
# Versions
CHIP_VER=2.0.2
R_VER=2.12.1
# Paths
EXEC_PATH=${PWD}
INST_PATH=/opt
CHIP_PATH=${INST_PATH}/chipster
TOOLS_PATH=${CHIP_PATH}/tools
TMPDIR_PATH=/tmp/install
# Misc
USERNAME=chipster
GROUPNAME=chipster

## Create tmpdir
rm -rf ${TMPDIR_PATH}/
mkdir ${TMPDIR_PATH}/

## Chipster:

# Install Chipster
cd ${TMPDIR_PATH}/
curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/versions/${CHIP_VER}/chipster-${CHIP_VER}.tar.gz | tar -xz
mv chipster/ ${CHIP_PATH}/

# Make some config "corrections"
sed -i'~' "s:/opt/chipster/:${CHIP_PATH}/:" ${CHIP_PATH}/comp/conf/environment.xml
sed -i'~' "s:/opt/chipster/:${CHIP_PATH}/:" ${CHIP_PATH}/comp/conf/runtimes.xml
sed -i'~' 's;http://www.bioconductor.org;http://mirrors.ebi.ac.uk/bioconductor/;' ${CHIP_PATH}/comp/conf/environment.xml
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

##############################################
# Install external applications and datasets #
##############################################

## In root:

# MACS, Artistic license
# part 1
cd ${TMPDIR_PATH}/
wget -nv http://liulab.dfci.harvard.edu/MACS/deb/macs_1.4.1.deb
dpkg -i macs_1.4.1.deb
rm macs_1.4.1.deb

# HTSeq, GPL v3 or later
# part 1
pip install HTSeq==0.5.3p3

## In tools:
 
if [ $mode == "devel" -a $build_tools == "yes" ]
then
  ## R:
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R-${R_VER}-vmbin/R-${R_VER}.tar.gz | tar -xz -C ${TOOLS_PATH}/  
  ln -s R-${R_VER} ${TOOLS_PATH}/R
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.12.1-vmbin/library/FruitFlyAgilent.db-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.12.1/lib64/R/library/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.12.1-vmbin/library/hgug4851a.db-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.12.1/lib64/R/library/
  # maSigPro is included in the R-2.12.1-vmbin.tar.gz
        
                        
#  cd ${TMPDIR_PATH}/
#  curl -s http://ftp.sunet.se/pub/lang/CRAN/src/base/R-2/R-${R_VER}.tar.gz | tar -xz
#  cd R-${R_VER}/
#  ## Fix for "/opt/chipster/tools/R-2.12.1/lib64/R/lib/libRlapack.so: undefined symbol: _gfortran_compare_string"
#  sed -i '/Rlapack_la_LIBADD =/ s/@DYLIB_UNDEFINED_ALLOWED_FALSE@//' src/modules/lapack/Makefile.in
#  export MAKEFLAGS=-j
#  #LIBnn=lib
#  ./configure --prefix=${TOOLS_PATH}/R-${R_VER}
#  make
#  make install
#  echo 'MAKEFLAGS=-j' > ${TOOLS_PATH}/R-${R_VER}/lib64/R/etc/Makevars.site # (could also be $HOME/.R/Makevars)
#  cd ../
#  rm -rf R-${R_VER}/
#  ln -s R-${R_VER} ${TOOLS_PATH}/R
#
#  ## R Libraries:
#  cd ${CHIP_PATH}/
#  echo | ./setup.sh

  # Add R package outside of setup.sh/environment.xml
  cd ${TMPDIR_PATH}/
  wget http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/nz131a520662fcdf.tar.gz
  ${TOOLS_PATH}/R/bin/R CMD INSTALL nz131a520662fcdf.tar.gz
  rm nz131a520662fcdf.tar.gz
  
  
      
  ## R-2.14:
  R_VER=2.14.1
  cd ${TMPDIR_PATH}/
  #wget -nv http://ftp.sunet.se/pub/lang/CRAN/src/base/R-2/R-${R_VER}.tar.gz
  #tar -xzf R-${R_VER}.tar.gz
  curl -s http://ftp.sunet.se/pub/lang/CRAN/src/base/R-2/R-${R_VER}.tar.gz | tar -xz
  cd R-${R_VER}/
  ## Fix for "/opt/chipster/tools/R-2.12.1/lib64/R/lib/libRlapack.so: undefined symbol: _gfortran_compare_string"
  #sed -i '/Rlapack_la_LIBADD =/ s/@DYLIB_UNDEFINED_ALLOWED_FALSE@//' src/modules/lapack/Makefile.in
  export MAKEFLAGS=-j
  #LIBnn=lib
  ./configure --prefix=${TOOLS_PATH}/R-${R_VER}
  make
  make install
  echo 'MAKEFLAGS=-j' > ${TOOLS_PATH}/R-${R_VER}/lib64/R/etc/Makevars.site # (could also be $HOME/.R/Makevars)
  cd ../
  rm -rf R-${R_VER}/
  ${TOOLS_PATH}/R-${R_VER}/bin/Rscript --vanilla ${CHIP_PATH}/comp/modules/admin/R-2.14/install-libs.R   


  ## R-2.15:
  R_VER=2.15.1
  cd ${TMPDIR_PATH}/
  curl -s http://ftp.sunet.se/pub/lang/CRAN/src/base/R-2/R-${R_VER}.tar.gz | tar -xz
  cd R-${R_VER}/
  export MAKEFLAGS=-j
  ./configure --prefix=${TOOLS_PATH}/R-${R_VER}
  make
  make install
  echo 'MAKEFLAGS=-j' > ${TOOLS_PATH}/R-${R_VER}/lib64/R/etc/Makevars.site # (could also be $HOME/.R/Makevars)
  cd ../
  rm -rf R-${R_VER}/
  ${TOOLS_PATH}/R-${R_VER}/bin/Rscript --vanilla ${CHIP_PATH}/comp/modules/admin/R-2.15/install-libs.R   

  # could also use the package from nic
  #curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-${R_VER}-vmbin/R-${R_VER}.tar.gz | tar -xz -C ${TOOLS_PATH}/  


  ## R-2.15.1_bioc-2.11	
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1_bioc-2.11/R-2.15.1_bioc-2.11-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/  
	

  ## External apps:

  # Link tool admin scripts from Chipster installation
  mkdir ${TOOLS_PATH}/admin/
  ln -s ${CHIP_PATH}/comp/modules/ngs/admin ${TOOLS_PATH}/admin/ngs
  
  # Weeder, custom license, according to developers VM bundling is ok
  cd ${TMPDIR_PATH}/
  curl -s http://159.149.109.9/modtools/downloads/weeder1.4.2.tar.gz | tar -xz
  cd Weeder1.4.2/
  ./compileall
  cd ../
  mkdir ${TOOLS_PATH}/weeder/
  mv Weeder1.4.2/ ${TOOLS_PATH}/weeder/

  # ClusterBuster, no license
  cd ${TMPDIR_PATH}/
  curl -s http://zlab.bu.edu/~mfrith/downloads/cbust-src.tar.gz | tar -xz
  cd cbust-src/
  make
  mkdir ${TOOLS_PATH}/ClusterBuster/
  mv cbust ${TOOLS_PATH}/ClusterBuster/
  cd ../
  rm -rf cbust-src/

  # Jaspar, no license
  cd ${TMPDIR_PATH}/
  wget -nv http://zlab.bu.edu/clover/jaspar2
  mv jaspar2 ${TOOLS_PATH}/ClusterBuster/jaspar2005core.txt

  # Promoter sequence files, license?
  cd ${TMPDIR_PATH}/
  mkdir ${TOOLS_PATH}/weeder/seqs/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/weeder_seqs/All_Weeder_sequence_files_v1.tar.gz | tar -xz -C ${TOOLS_PATH}/weeder/seqs/

  # BEDTools, GNU GPL v2
  cd ${TMPDIR_PATH}/
  curl -s http://bedtools.googlecode.com/files/BEDTools.v2.12.0.tar.gz | tar -xz
  cd BEDTools-Version-2.12.0
  make clean
  make all
  cd ../
  mv BEDTools-Version-2.12.0/ ${TOOLS_PATH}/
  ln -s BEDTools-Version-2.12.0 ${TOOLS_PATH}/bedtools

  # MACS, Artistic license
  # part 2
  cd ${TMPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/macs/
  ln -s /usr/bin/macs14 ${TOOLS_PATH}/macs/macs14

  # SAM tools, BSD License, MIT License
  cd ${TMPDIR_PATH}/
  #curl -sL http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2/download | tar -xj
  #cd samtools-0.1.18/
  #make
  #cd ../
  #mv samtools-0.1.18/ ${TOOLS_PATH}
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/samtools-0.1.18.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s samtools-0.1.18 ${TOOLS_PATH}/samtools

  # Bowtie, Artistic License
  cd ${TMPDIR_PATH}/
  wget -nv -O bowtie-0.12.7-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download
  unzip -q bowtie-0.12.7-linux-x86_64.zip
  mv bowtie-0.12.7/ ${TOOLS_PATH}
  ln -s bowtie-0.12.7 ${TOOLS_PATH}/bowtie
  rm bowtie-0.12.7-linux-x86_64.zip

  # Bowtie 2, Artistic License
  cd ${TMPDIR_PATH}/
  wget -nv -O bowtie2-2.0.0-beta7-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.0.0-beta7/bowtie2-2.0.0-beta7-linux-x86_64.zip/download
  unzip -q bowtie2-2.0.0-beta7-linux-x86_64.zip
  mv bowtie2-2.0.0-beta7 ${TOOLS_PATH}
  ln -s bowtie2-2.0.0-beta7 ${TOOLS_PATH}/bowtie2
  rm bowtie2-2.0.0-beta7-linux-x86_64.zip

	# Fasta files
  cd ${TMPDIR_PATH}/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_Halorubrum_lacusprofundi_ATCC_49239.tar.gz | tar -xz -C ${TOOLS_PATH}/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_e_coli.tar.gz | tar -xz -C ${TOOLS_PATH}/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_hg19.tar.gz | tar -xz -C ${TOOLS_PATH}/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_mm10.tar.gz | tar -xz -C ${TOOLS_PATH}/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_mm9.tar.gz | tar -xz -C ${TOOLS_PATH}/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_ovis_aries_texel.tar.gz | tar -xz -C ${TOOLS_PATH}/

	# Fasta files, nochr
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_hg19.tar.gz | tar -xz -C ${TOOLS_PATH}/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_mm10.tar.gz | tar -xz -C ${TOOLS_PATH}/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_mm9.tar.gz | tar -xz -C ${TOOLS_PATH}/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/

  # Bowtie indexes, built for Chipster
  cd ${TMPDIR_PATH}/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_Gasterosteus_aculeatus.BROADS1.67.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_Halorubrum_lacusprofundi_ATCC_49239.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_athaliana.TAIR10.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_canFam2.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_e_coli.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_hg19.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_miRBase18_mmu_matureT.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_mm10.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_mm9.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_ovis_aries_texel.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	ln -s -t ${TOOLS_PATH}/bowtie/indexes ../../genomes/fasta/nochr

	# Bowtie2 indexes, built for Chipster
  cd ${TMPDIR_PATH}/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_athaliana.TAIR10.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_canFam2.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_e_coli.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Halorubrum_lacusprofundi_ATCC_49239.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_hg19.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_mm9.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_mm10.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_ovis_aries_texel.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
	ln -s -t ${TOOLS_PATH}/bowtie2/indexes ../../genomes/fasta/nochr
				
  # FastQC, GPL v3 or later
  cd ${TMPDIR_PATH}/
  wget -nv http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/fastqc_v0.10.0.zip
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

  # Cufflinks, Boost License
  cd ${TMPDIR_PATH}/
  curl -s http://cufflinks.cbcb.umd.edu/downloads/cufflinks-1.0.3.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/  
  ln -s cufflinks-1.0.3.Linux_x86_64 ${TOOLS_PATH}/cufflinks
  curl -s http://cufflinks.cbcb.umd.edu/downloads/cufflinks-2.0.2.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/
 
  # Tophat, The Artistic License
  cd ${TMPDIR_PATH}/
  curl -s http://tophat.cbcb.umd.edu/downloads/tophat-1.3.2.Linux_x86_64.tar.gz | tar -xz
  mv tophat-1.3.2.Linux_x86_64 ${TOOLS_PATH}/
  ln -s tophat-1.3.2.Linux_x86_64 ${TOOLS_PATH}/tophat

  # Tophat 2, The Artistic License
  cd ${TMPDIR_PATH}/
  curl -s http://tophat.cbcb.umd.edu/downloads/tophat-2.0.4.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s tophat-2.0.4.Linux_x86_64 ${TOOLS_PATH}/tophat2

  # BWA, GPL v3 or later, MIT License
  cd ${TMPDIR_PATH}/
  curl -sL http://sourceforge.net/projects/bio-bwa/files/bwa-0.6.1.tar.bz2/download | tar -xj
  cd bwa-0.6.1/
  make
  cd ../
  mv bwa-0.6.1/ ${TOOLS_PATH}/
  ln -s bwa-0.6.1 ${TOOLS_PATH}/bwa

  # BWA index check
  curl -sL http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/misc/check_bwa_index.sh > ${TOOLS_PATH}/bwa/check_bwa_index.sh
  chmod 755 ${TOOLS_PATH}/bwa/check_bwa_index.sh

  # Fastx links
  mkdir -p ${TOOLS_PATH}/fastx/bin/
  ln -s /usr/bin/fasta_* ${TOOLS_PATH}/fastx/bin/
  ln -s /usr/bin/fastq_* ${TOOLS_PATH}/fastx/bin/
  ln -s /usr/bin/fastx_* ${TOOLS_PATH}/fastx/bin/

  # CanGEM probe mapping data
  cd ${TMPDIR_PATH}/
  mkdir ${TOOLS_PATH}/CanGEM/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/CanGEM_probe_mappings/All_CanGEM_probe_mappings_v1.tar.gz | tar -xz -C ${TOOLS_PATH}/CanGEM/

  # Genome data for tools
  cd ${TMPDIR_PATH}/
  mkdir ${TOOLS_PATH}/genomes/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes_for_tools/All_genomes_for_tools_v1.tar.gz | tar -xz -C ${TOOLS_PATH}/genomes/

  # GTF gene data for tools
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Homo_sapiens.GRCh37.68.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Mus_musculus.GRCm38.68.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Rattus_norvegicus.RGSC3.4.68.tar.gz | tar -xz -C ${TOOLS_PATH}/

  # miRNA mapping data
  cd ${TMPDIR_PATH}/
  mkdir ${TOOLS_PATH}/miRNA_mappings/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/miRNA_mappings/All_miRNA_mappings_v1.tar.gz | tar -xz -C ${TOOLS_PATH}/miRNA_mappings/

  # Genome variant databases
  cd ${TMPDIR_PATH}/
  mkdir ${TOOLS_PATH}/DGV/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomic_variant_dbs/All_genomic_variant_dbs_v1.tar.gz | tar -xz -C ${TOOLS_PATH}/DGV/

  # bwa indexes, built for Chipster
  cd ${TMPDIR_PATH}/
  mkdir ${TOOLS_PATH}/bwa_indexes/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_hg9.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_mm9.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_mmu_miRB17mature.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_ovis_aries_texel.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_mm10.tar.gz | tar -xz -C ${TOOLS_PATH}/


  # Data for CNA-seq tools (produced by Ilari Scheinin)
  cd ${TMPDIR_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/CNA_seq/MPScall.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/CNA_seq/FREEC_Linux64.tar.gz | tar -xz -C ${TOOLS_PATH}/

  # prinseq
  cd ${TMPDIR_PATH}/
  curl -sL http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.19.3.tar.gz/download | tar -xz
  chmod a+x prinseq-lite-0.19.3/prinseq-lite.pl
  chmod a+x prinseq-lite-0.19.3/prinseq-graphs.pl
  mv prinseq-lite-0.19.3 ${TOOLS_PATH}/
  ln -s prinseq-lite-0.19.3 ${TOOLS_PATH}/prinseq

  # Genome annotations for genome browser
  cd ${TMPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/All_genomes_for_browser_v2.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/genomebrowser_fasta_Ovis_aries.Oar_v3.1.dna.toplevel.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  wget -O ${TOOLS_PATH}/genomebrowser/annotations/contents2.txt http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/contents2.txt
  
  # DEXSeq
	cd ${TMPDIR_PATH}/
	#	curl -sL http://www.bioconductor.org/packages/release/bioc/src/contrib/DEXSeq_1.2.1.tar.gz | tar -xz
	curl -sL http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/DEXSeq_1.2.1.tar.gz | tar -xz
	mkdir ${TOOLS_PATH}/dexseq-exoncounts
	cp DEXSeq/inst/python_scripts/dexseq_count.py ${TOOLS_PATH}/dexseq-exoncounts  
	cp DEXSeq/inst/python_scripts/dexseq_prepare_annotation.py ${TOOLS_PATH}/dexseq-exoncounts
	rm -rf DEXSeq	

 	# vcftools, GPLv3
    cd ${TMPDIR_PATH}/
    curl -sL http://sourceforge.net/projects/vcftools/files/vcftools_0.1.9.tar.gz/download| tar -xz
    cd vcftools_0.1.9/
    make
    cd ../
    mv vcftools_0.1.9/ ${TOOLS_PATH}/
    ln -s vcftools_0.1.9 ${TOOLS_PATH}/vcftools

  # GATK, MIT
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/GenomeAnalysisTKLite-latest.tar.bz2 | tar -xj -C ${TOOLS_PATH}/
  ln -s GenomeAnalysisTKLite-2.1-11-gfb37f33 ${TOOLS_PATH}/GenomeAnalysisTK2

 	 	 	 	  	 	 	 	  	 	 	 	 
  ## Create checksums
  cd ${TOOLS_PATH}/
  find . '!' -type d '!' -type l -print0 | xargs -0 sha256sum >> tools.sha256sum
fi

## Clean up:
chown -R -h ${USERNAME}:${GROUPNAME} ${CHIP_PATH}/
chown -R ${USERNAME}:${GROUPNAME} /mnt/tools/
ls -lah ${TMPDIR_PATH}/ # Debug, list uncleaned mess
rm -rf ${TMPDIR_PATH}/

## Configure Chipster:
# Configure chipster ports and hostname

## Init.d:
ln -s ${CHIP_PATH}/chipster /etc/init.d/chipster
update-rc.d chipster defaults 30 70 # start 99 3 5 . stop 1 0 1 2 4 6 .

ln -s ${CHIP_PATH}/activemq/bin/linux-x86-64/activemq /etc/init.d/chipster-activemq
ln -s ${CHIP_PATH}/comp/bin/linux-x86-64/chipster-comp /etc/init.d/chipster-comp
ln -s ${CHIP_PATH}/auth/bin/linux-x86-64/chipster-auth /etc/init.d/chipster-auth
ln -s ${CHIP_PATH}/fileserver/bin/linux-x86-64/chipster-fileserver /etc/init.d/chipster-fileserver
ln -s ${CHIP_PATH}/webstart/bin/linux-x86-64/chipster-webstart /etc/init.d/chipster-webstart
ln -s ${CHIP_PATH}/manager/bin/linux-x86-64/chipster-manager /etc/init.d/chipster-manager
#update-rc.d chipster-activemq defaults
#update-rc.d chipster-comp defaults
#update-rc.d chipster-auth defaults
#update-rc.d chipster-fileserver defaults
#update-rc.d chipster-webstart defaults
#update-rc.d chipster-manager defaults
