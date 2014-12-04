#!/bin/bash

# This script updates to latest version. Updates between minor versions should be smooth and
# automatic, where as updates between major versions can require some manual steps afterwards
# if some specific local customizations were in place.
# This update mechanism has been available since 2.0.2.

# Latest version, matching tar-packages must be available 
##
LATEST_VERSION=2.12.1
R_VERSION=3.0.2

# Exit immediately if some command fails
set -e

# Helper functions
function compare_versions()
{
  VERSION1_ARR=( `echo $1 | tr "." "\n"` )
  VERSION2_ARR=( `echo $2 | tr "." "\n"`)

  # Check main version  
  if [ ${VERSION1_ARR[0]} -lt ${VERSION2_ARR[0]} ] ; then
    COMPARE_VERSIONS_RESULT=1
    return
  fi
  if [ ${VERSION1_ARR[0]} -gt ${VERSION2_ARR[0]} ] ; then
    COMPARE_VERSIONS_RESULT=-1
    return
  fi

  # Check major version  
  if [ ${VERSION1_ARR[1]} -lt ${VERSION2_ARR[1]} ] ; then
    COMPARE_VERSIONS_RESULT=1
    return
  fi
  if [ ${VERSION1_ARR[1]} -gt ${VERSION2_ARR[1]} ] ; then
    COMPARE_VERSIONS_RESULT=-1
    return
  fi
  
  # Check minor version
  if [ ${VERSION1_ARR[2]} -lt ${VERSION2_ARR[2]} ] ; then
    COMPARE_VERSIONS_RESULT=1
    return
  fi
  if [ ${VERSION1_ARR[2]} -gt ${VERSION2_ARR[2]} ] ; then
    COMPARE_VERSIONS_RESULT=-1
    return
  fi
  
  COMPARE_VERSIONS_RESULT=0
}


function compare_to_current()
{
    compare_versions $1 $CURRENT_VERSION
    CURRENT_COMPARED=$COMPARE_VERSIONS_RESULT
}

function compare_to_latest()
{
    compare_versions $1 $LATEST_VERSION
    LATEST_COMPARED=$COMPARE_VERSIONS_RESULT
}

function compare_to_current_and_latest()
{
    compare_to_current $1
    compare_to_latest $1
}


function info_older_than_252()
{
  echo "IMPORTANT!"
  echo "The latest version of Chipster is $LATEST_VERSION and your current version is $CURRENT_VERSION."
  echo "The update to Chipster $LATEST_VERSION needs to be done in two phases: first update to 2.5.2"
  echo "and after that update to $LATEST_VERSION."
  echo "The update will next proceed to update Chipster to 2.5.2. After this update is finished,"
  echo "you need to run ./update.sh again to continue to update to Chipster $LATEST_VERSION."
}

function info_before_260()
{
  echo "IMPORTANT! Please read this!"
  echo ""
  echo "There are two main disk images in the Chipster virtual machine: root.vmdk and tools.vmdk."
  echo "Root.vmdk contains operating system and Chipster binaries"
  echo "and configuration files. The size of the root.vmdk is less than 2 GB at the moment."
  echo "Tools.vmdk contains analysis tools binaries, libraries and"
  echo "annotations. The size of the tools.vmdk is a bit less than 200 GB at the moment."
  echo ""
  echo "In the original Chipster VM the maximum size of the tools.vmdk was about 188 GB."
  echo "This maximum was later increased to 512 GB."
  echo ""
  echo "Since Chipster 2.6.0, the original tools image has become too small to contain all"
  echo "the necessary files for Chipster. Therefore, depending on the time you first downloaded"
  echo "Chipster VM, your tools image may be too small for Chipster 2.6".
  echo ""
  echo "To find out whether this is the case, check the 'Size' column in the next"
  echo "few lines."
  echo ""
  df -h /opt/chipster/tools
  echo ""
  echo "If the size is 512 GB, you can continue to update to 2.6. If the size is 188 GB, you"
  echo "need to download Chipster 2.5.2 tools.vmdk and replace your current tools.vmdk with it."
  echo "After that, you can continue to update to Chipster 2.6 and later."
  echo ""
  echo "You can also always download the whole Chipster VM, but by downloading only tools.vmdk, you"
  echo "can keep all the configurations which are located at root.vmdk."
  echo ""
}


# Note for Chipster 3.x
echo ""
echo "IMPORTANT!"
echo ""
echo "You are running Chipster 2.x while the latest Chipster is 3.x."
echo ""
echo "With this update script you can only update Chipster to the latest 2.x version."
echo ""
echo "It is highly recommended to get the Chipster 3.0, by downloading the new virtual machine from"
echo "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/virtual_machines/"
echo ""
echo "Continue with the update?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) break;;
        No )  exit;;
    esac
done
echo ""



# Make sure user has sudo rights
echo ""
echo "Some parts of the update may need root privileges. These parts are run using sudo."
echo "Testing permission to use sudo..."
if [ "$(sudo whoami)" != 'root' ]; then echo 'You need sudo rights to run the update script, aborting.'; 
exit 1; fi
echo "Sudo ok"
echo ""

# Detect current version
CURRENT_VERSION=`ls -1 shared/lib | grep ^chipster-[0-9\\.]*.jar | gawk 'match($0, "chipster-([0-9\\\\.]*).jar", g) {print g[1]}'`


## before 2.6 specials                                                                                                                                                                                               
compare_to_current "2.5.2"

# current is older than 2.5.2                                                                                                                                                                               
if [ $CURRENT_COMPARED -lt 0 ]
then
    info_older_than_252
    END_MESSAGE="IMPORTANT: Run ./update.sh again to update to Chipster $LATEST_VERSION"
    LATEST_VERSION=2.5.2
    echo "Continue with the update?"
    select yn in "Yes" "No"; do
        case $yn in
            Yes ) break;;
            No )  exit;;
        esac
    done
    echo ""

# current is 2.5.2                                                                                                                                                                                          
elif [ $CURRENT_COMPARED -eq 0 ]
then
  info_before_260
  SELECT_UPDATE_TO_LATEST="I have the 512 GB tools image and want to continue the update now"
  SELECT_UPDATE_LATER="I have the old 188 GB tools image, so I'll exit now, download the Chipster 2.5.2 tools.vmdk and then run update.sh again"

  options=("${SELECT_UPDATE_TO_LATEST}" "${SELECT_UPDATE_LATER}")
  PS3='Please enter your choice: '
  select opt in "${options[@]}"
  do
    case $opt in
        $SELECT_UPDATE_TO_LATEST)
            echo "Continuing update normally"
            break
            ;;
        $SELECT_UPDATE_LATER)
            echo "Quitting without updating"
            echo ""
            echo "Next steps:"
            echo "- Download and Chipster 2.5.2 tools.vmdk http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/virtual_machines/2.5.2/VMware/tools.vmdk"
            echo "- Shutdown Chipster virtual machine"
            echo "- Replace your current tools.vmdk with the downloaded tools.vmdk"
            echo "- Start Chipster virtual machine"
            echo "- Run update.sh again"
            echo ""
            exit;
            ;;
        *) echo invalid option;;
    esac
  done
fi





## normal 2.6 and later update

# Check current version
echo Detected version $CURRENT_VERSION
compare_to_current "$LATEST_VERSION"
if [ $CURRENT_COMPARED -gt 0 ] ; then 
  echo "Update error: current version $CURRENT_VERSION is newer than latest $LATEST_VERSION"
  exit 1
fi
if [ $CURRENT_COMPARED -eq 0 ] ; then 
  echo "Already at latest version, nothing needs to be updated"
  exit
fi
echo "Will update to version $LATEST_VERSION"
echo ""

# EMBOSS warning for 2.12
echo "IMPORTANT!"
echo ""
echo "When updating to Chipster 2.12, it is highly recommended to download the full virtual machine"
echo "images from http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/virtual_machines/2.12.0/"
echo ""
echo "If you update to 2.12 using this script, some of the EMBOSS tools may not work correctly."
echo ""
echo "Continue with the update?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) break;;
        No )  exit;;
    esac
done
echo ""


# Confirm update
echo "Update will start next. It can take several hours, depending on your network connection"
echo "IMPORTANT: Stop the Chipster service before proceeding with the update: 'service chipster stop'"
echo "Do you wish to proceed with the update?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) echo "** Update started"; break;;
        No ) echo "** Update aborted"; exit;;
    esac
done

# Start update
INST_PATH=/opt
CHIP_PATH=${INST_PATH}/chipster
TOOLS_PATH=${CHIP_PATH}/tools
TMPDIR_PATH=/tmp/chipster-install-temp

# Create temp dir
rm -rf ${TMPDIR_PATH}
mkdir ${TMPDIR_PATH}

# Create backup dir
TIMESTAMP=$(date +"%Y%m%d%H%M")
BACKUPDIR_PATH=/tmp/$TIMESTAMP-$RANDOM
while [ -d "$BACKUPDIR_PATH" ] ; do
  BACKUPDIR_PATH=/tmp/$TIMESTAMP-$RANDOM
done
mkdir ${BACKUPDIR_PATH}

#######################################
# VERSION SPECIFIC ENTRIES START HERE #
# (ADD NEW ENTRIES TO THE END)        #
#######################################

# 2.0.3
compare_to_current_and_latest "2.0.3"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Installing mm10 bowtie indexes"
  curl -L http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_mm10.tar.gz  | tar -xz -C ${TOOLS_PATH}/bowtie/
fi

# 2.1.0 
compare_to_current_and_latest "2.1.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

    echo "** Updating prinseq"
    cd ${TMPDIR_PATH}/
    curl -L http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.19.3.tar.gz/download | tar -xz
    chmod a+x prinseq-lite-0.19.3/prinseq-lite.pl
    chmod a+x prinseq-lite-0.19.3/prinseq-graphs.pl
    rm -rf ${TOOLS_PATH}/prinseq*
    mv prinseq-lite-0.19.3 ${TOOLS_PATH}/
    ln -s prinseq-lite-0.19.3 ${TOOLS_PATH}/prinseq

    echo "Adding nz131a520662fcdf"
    cd ${TMPDIR_PATH}/
    wget http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/nz131a520662fcdf.tar.gz
    ${TOOLS_PATH}/R/bin/R CMD INSTALL nz131a520662fcdf.tar.gz
    rm nz131a520662fcdf.tar.gz
    
    echo "Installing vcftools"
    cd ${TMPDIR_PATH}/
    curl -L http://sourceforge.net/projects/vcftools/files/vcftools_0.1.9.tar.gz/download | tar -xz
    cd vcftools_0.1.9/
    make
    cd ../
    rm -rf ${TOOLS_PATH}/vcftools*
    mv vcftools_0.1.9/ ${TOOLS_PATH}/
    ln -s vcftools_0.1.9 ${TOOLS_PATH}/vcftools
    
    echo "Updating genome browser annotations"
    mv -b ${TOOLS_PATH}/genomebrowser/annotations ${BACKUPDIR_PATH}/
    mkdir ${TOOLS_PATH}/genomebrowser/annotations # not typically needed, but the tar package is a bit stupid in this case
    curl -L http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/All_genomes_for_browser_v2.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  
    echo "Installing R-2.15"
    curl -L http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1-vmbin/R-2.15.1.tar.gz | tar -xz -C ${TOOLS_PATH}/
        
fi


# 2.2.0 
compare_to_current_and_latest "2.2.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Updating samtools"
  cd ${TMPDIR_PATH}/
  mv -b ${TOOLS_PATH}/samtools-0.1.13/ ${BACKUPDIR_PATH}/
  rm ${TOOLS_PATH}/samtools
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/samtools-0.1.18.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s samtools-0.1.18 ${TOOLS_PATH}/samtools

  echo "** Relocating fasta files"
  mkdir ${TOOLS_PATH}/genomes/fasta
  rm -f ${TOOLS_PATH}/bowtie/indexes/nochr/*.fai
  mv ${TOOLS_PATH}/bowtie/indexes/*.fa ${TOOLS_PATH}/genomes/fasta/
  mv ${TOOLS_PATH}/bowtie/indexes/nochr ${TOOLS_PATH}/genomes/fasta/
  ln -s -t ${TOOLS_PATH}/bowtie/indexes ../../genomes/fasta/e_coli.fa 
  ln -s -t ${TOOLS_PATH}/bowtie/indexes ../../genomes/fasta/Halorubrum_lacusprofundi_ATCC_49239.fa 
  ln -s -t ${TOOLS_PATH}/bowtie/indexes ../../genomes/fasta/hg19.fa 
  ln -s -t ${TOOLS_PATH}/bowtie/indexes ../../genomes/fasta/mm10.fa
  ln -s -t ${TOOLS_PATH}/bowtie/indexes ../../genomes/fasta/mm9.fa
  ln -s -t ${TOOLS_PATH}/bowtie/indexes ../../genomes/fasta/Phytophthora_infestans1_1.12.fa
  ln -s -t ${TOOLS_PATH}/bowtie/indexes ../../genomes/fasta/Populus_trichocarpa.JGI2.0.12.fa
  ln -s -t ${TOOLS_PATH}/bowtie/indexes ../../genomes/fasta/rn4.fa
  ln -s -t ${TOOLS_PATH}/bowtie/indexes ../../genomes/fasta/nochr
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_mm10.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_mm10.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Installing bowtie2"
  cd ${TMPDIR_PATH}/
  wget -nv -O bowtie2-2.0.0-beta7-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.0.0-beta7/bowtie2-2.0.0-beta7-linux-x86_64.zip/download
  unzip -q bowtie2-2.0.0-beta7-linux-x86_64.zip
  mv bowtie2-2.0.0-beta7 ${TOOLS_PATH}
  ln -s bowtie2-2.0.0-beta7 ${TOOLS_PATH}/bowtie2
  rm bowtie2-2.0.0-beta7-linux-x86_64.zip

  echo "** Installing bowtie2 indexes"
  cd ${TMPDIR_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_athaliana.TAIR10.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_canFam2.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_e_coli.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Halorubrum_lacusprofundi_ATCC_49239.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_hg19.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_mm10.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Phytophthora_infestans1_1.12.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Populus_trichocarpa.JGI2.0.12.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_rn4.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_saprolegnia_parasitica_cbs_223.65_2_contigs.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  ln -s -t ${TOOLS_PATH}/bowtie2/indexes ../../genomes/fasta/nochr
  
  echo "** Installing tophat2"
  cd ${TMPDIR_PATH}/
  curl -s http://tophat.cbcb.umd.edu/downloads/tophat-2.0.4.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s tophat-2.0.4.Linux_x86_64 ${TOOLS_PATH}/tophat2

  echo "** Installing DEXSeq"
  cd ${TMPDIR_PATH}/
  curl -sL http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/DEXSeq_1.2.1.tar.gz | tar -xz
  mkdir ${TOOLS_PATH}/dexseq-exoncounts
  cp DEXSeq/inst/python_scripts/dexseq_count.py ${TOOLS_PATH}/dexseq-exoncounts  
  cp DEXSeq/inst/python_scripts/dexseq_prepare_annotation.py ${TOOLS_PATH}/dexseq-exoncounts
  rm -rf DEXSeq 

  echo "** Updating gtf files"
  mkdir -p ${BACKUPDIR_PATH}/genomes
  mv -b ${TOOLS_PATH}/genomes/Homo_sapiens.GRCh37.62.chr.gtf ${BACKUPDIR_PATH}/
  mv -b ${TOOLS_PATH}/genomes/Mus_musculus.NCBIM37.62.chr.gtf ${BACKUPDIR_PATH}/
  mv -b ${TOOLS_PATH}/genomes/Rattus_norvegicus.RGSC3.4.62.chr.gtf ${BACKUPDIR_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Homo_sapiens.GRCh37.68.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Mus_musculus.GRCm38.68.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Rattus_norvegicus.RGSC3.4.68.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Installing cufflinks2" 
  cd ${TMPDIR_PATH}/
  curl -s http://cufflinks.cbcb.umd.edu/downloads/cufflinks-2.0.2.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Fixing vcftools"
  cd ${TMPDIR_PATH}/
  rm -rf ${TOOLS_PATH}/vcftools*
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/vcftools_0.1.9.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s vcftools_0.1.9 ${TOOLS_PATH}/vcftools
                        
fi

# 2.2.1
compare_to_current_and_latest "2.2.1"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 
  echo "** Tools image is uptodate"                      
fi

# 2.2.2
compare_to_current_and_latest "2.2.2"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 
  echo "** Removing obsolete files"                      
  rm -f ${TOOLS_PATH}/admin/ngs
  rm -rf ${TOOLS_PATH}/tophat-1.3.0.Linux_x86_64  
  rm -f ${TOOLS_PATH}/bwa_indexes/*.fa
  
  echo "** Adding sheep genome"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/genomebrowser_fasta_Ovis_aries.Oar_v3.1.dna.toplevel.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  wget -O ${TOOLS_PATH}/genomebrowser/annotations/contents2.txt http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/contents2.txt
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_ovis_aries_texel.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_ovis_aries_texel.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_ovis_aries_texel.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_ovis_aries_texel.tar.gz | tar -xz -C ${TOOLS_PATH}/
  
  echo "** Installing R library maSigPro"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.12.1-vmbin/library/maSigPro-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.12.1/lib64/R/library/
  
  echo "** Updating R-2.15"
  mv -b ${TOOLS_PATH}/R-2.15.1 ${BACKUPDIR_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1-vmbin/R-2.15.1.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Installing GATK"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/GenomeAnalysisTKLite-latest.tar.bz2 | tar -xj -C ${TOOLS_PATH}/
  ln -s GenomeAnalysisTKLite-2.1-11-gfb37f33 ${TOOLS_PATH}/GenomeAnalysisTK2

fi

# 2.2.3
compare_to_current_and_latest "2.2.3"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Installing R library FruitFlyAgilent.db"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.12.1-vmbin/library/FruitFlyAgilent.db-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.12.1/lib64/R/library/

  echo "** Removing obsolete genomes and indexes"
  rm -f ${TOOLS_PATH}/genomes/fasta/Phytophthora_infestans1_1.12.fa
  rm -f ${TOOLS_PATH}/genomes/fasta/Populus_trichocarpa.JGI2.0.12.fa
  rm -f ${TOOLS_PATH}/bowtie/indexes/Phytophthora_infestans1_1.12.*
  rm -f ${TOOLS_PATH}/bowtie/indexes/Populus_trichocarpa.JGI2.0.12.*
  rm -f ${TOOLS_PATH}/bowtie/indexes/saprolegnia_parasitica_cbs_223.65_2_contigs.*
  rm -f ${TOOLS_PATH}/bowtie2/indexes/Phytophthora_infestans1_1.12.*
  rm -f ${TOOLS_PATH}/bowtie2/indexes/Populus_trichocarpa.JGI2.0.12.*
  rm -f ${TOOLS_PATH}/bowtie2/indexes/saprolegnia_parasitica_cbs_223.65_2_contigs.*

fi


# 2.2.4
compare_to_current_and_latest "2.2.4"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Installing mm10 bwa index"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_mm10.tar.gz | tar -xz -C ${TOOLS_PATH}/
  
  echo "** Updating R library FruitFlyAgilent.db"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.12.1-vmbin/library/FruitFlyAgilent.db-vmbin.tar.gz | tar -xz --overwrite -C ${TOOLS_PATH}/R-2.12.1/lib64/R/library/

  echo "** Installing R library hgug4851a.db"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.12.1-vmbin/library/hgug4851a.db-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.12.1/lib64/R/library/
        
fi


# 2.3.0
compare_to_current_and_latest "2.3.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Installing R-2.15 with Bioconductor 2.11"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1_bioc-2.11/R-2.15.1_bioc-2.11-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Updating genome browser annotations"
  mv -b ${TOOLS_PATH}/genomebrowser/annotations ${BACKUPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/genomebrowser/annotations

  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Arabidopsis_lyrata.v.1.0.16.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Arabidopsis_thaliana.TAIR10.16.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Canis_familiaris.BROADD2.67.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Canis_familiaris.CanFam3.1.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Gasterosteus_aculeatus.BROADS1.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Homo_sapiens.GRCh37.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Homo_sapiens.NCBI36.54.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Mus_musculus.GRCm38.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Mus_musculus.NCBIM37.67.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Ovis_aries_v3.1.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Rattus_norvegicus.RGSC3.4.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Sus_scrofa.Sscrofa10.2.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Vitis_vinifera.IGGP_12x.16.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/Yersinia.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  wget -O ${TOOLS_PATH}/genomebrowser/annotations/contents2.txt http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.3.0/contents2.txt

fi

# 2.3.1
compare_to_current_and_latest "2.3.1"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Updating R library VariantAnnotation"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1_bioc-2.11/library/VariantAnnotation-vmbin.tar.gz | tar -xz --overwrite -C ${TOOLS_PATH}/R-2.15.1_bioc-2.11/lib64/R/library/

fi

# 2.3.2
compare_to_current_and_latest "2.3.2"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then

  echo "** Updating R-2.15 with Bioconductor 2.11"
  mv -b ${TOOLS_PATH}/R-2.15.1_bioc-2.11 ${BACKUPDIR_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1_bioc-2.11/R-2.15.1_bioc-2.11-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Installing Arabidopsis lyrata fasta, bwa index and bowtie2 index"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_Arabidopsis_lyrata.v.1.0.16.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_Arabidopsis_lyrata.v.1.0.16.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Arabidopsis_lyrata.v.1.0.16.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/

  echo "** Installing Sus scrofa fasta, bwa index and bowtie2 index"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_Sus_scrofa.Sscrofa10.2.69.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Sus_scrofa.Sscrofa10.2.69.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_Sus_scrofa.Sscrofa10.2.69.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Sus_scrofa.Sscrofa10.2.69.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/

  echo "** Creating link for cufflinks2"
  ln -s cufflinks-2.0.2.Linux_x86_64 ${TOOLS_PATH}/cufflinks2

fi

# 2.4.0
compare_to_current_and_latest "2.4.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Installing genome browser fasta files"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_N916Ysi.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_R1-RT.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Arabidopsis_thaliana.TAIR10.16.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Bos_taurus.UMD3.1.69.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Canis_familiaris.BROADD2.67.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Canis_familiaris.CanFam3.1.69.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Gallus_gallus.Gallus_gallus-4.0.pre.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Gasterosteus_aculeatus.BROADS1.69.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Homo_sapiens.NCBI36.54.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Human-MT.NC_012920.1.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Vitis_vinifera.IGGP_12x.16.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Updating genome browser annotations"
  mv -b ${TOOLS_PATH}/genomebrowser/annotations ${BACKUPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/genomebrowser/annotations

  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Arabidopsis_lyrata.v.1.0.16.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Arabidopsis_thaliana.TAIR10.16.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Bos_taurus.UMD3.1.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Canis_familiaris.BROADD2.67.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Canis_familiaris.CanFam3.1.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Gallus_gallus.Gallus_gallus-4.0.pre.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Gasterosteus_aculeatus.BROADS1.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Homo_sapiens.GRCh37.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Homo_sapiens.NCBI36.54.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Human-MT.NC_012920.1.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Mus_musculus.GRCm38.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Mus_musculus.NCBIM37.67.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/N916Ysi.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Ovis_aries_v3.1.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/R1-RT.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Rattus_norvegicus.RGSC3.4.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Sus_scrofa.Sscrofa10.2.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/Vitis_vinifera.IGGP_12x.16.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/

  wget -O ${TOOLS_PATH}/genomebrowser/annotations/contents2.txt http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.4.0/contents2.txt

  echo "** Updating activemq.xml"
  cp ${CHIP_PATH}/activemq/conf/activemq.xml ${BACKUPDIR_PATH}/
  sed -i'~' "s:<authorizationEntry admin=\"all\" read=\"managers\" topic=\"feedback-topic\" write=\"clients\"/>:<authorizationEntry admin=\"all\" read=\"authenticators\" topic=\"feedback-topic\" write=\"clients\"/>\n\t\t\t\t\t<authorizationEntry admin=\"all\" read=\"managers\" topic=\"authorised-feedback-topic\" write=\"authenticators\"/>:" ${CHIP_PATH}/activemq/conf/activemq.xml

fi


# 2.5.0
compare_to_current_and_latest "2.5.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Updating R-2.15, Bioconductor 2.11"
  mv -b ${TOOLS_PATH}/R-2.15.1_bioc-2.11 ${BACKUPDIR_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1_bioc-2.11/R-2.15.1_bioc-2.11-vmbin_v2.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Updating fasta files"
  rm -f ${TOOLS_PATH}/genomes/fasta/Arabidopsis_thaliana.TAIR10.16.dna.toplevel.fa
  rm -f ${TOOLS_PATH}/genomes/fasta/Bos_taurus.UMD3.1.69.dna.toplevel.fa
  rm -f ${TOOLS_PATH}/genomes/fasta/Canis_familiaris.BROADD2.67.dna.toplevel.fa
  rm -f ${TOOLS_PATH}/genomes/fasta/Canis_familiaris.CanFam3.1.69.dna.toplevel.fa
  rm -f ${TOOLS_PATH}/genomes/fasta/Gallus_gallus.Gallus_gallus-4.0.pre.fa
  rm -f ${TOOLS_PATH}/genomes/fasta/Gasterosteus_aculeatus.BROADS1.69.dna.toplevel.fa
  rm -f ${TOOLS_PATH}/genomes/fasta/Homo_sapiens.NCBI36.54.dna.toplevel.fa
  rm -f ${TOOLS_PATH}/genomes/fasta/Human-MT.NC_012920.1.fa
  rm -f ${TOOLS_PATH}/genomes/fasta/Vitis_vinifera.IGGP_12x.16.dna.toplevel.fa

  # New fasta packages with 'nochr' folder
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Arabidopsis_thaliana.TAIR10.17.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Bos_taurus.UMD3.1.70.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Canis_familiaris.BROADD2.67-v2.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Canis_familiaris.CanFam3.1.70.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Gallus_gallus.Gallus_gallus-4.0.pre-v2.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Gasterosteus_aculeatus.BROADS1.70.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Homo_sapiens.NCBI36.54-v2.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Human-MT.NC_012920.1-v2.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Vitis_vinifera.IGGP_12x.17.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Updating genome browser annotations"
  mv -b ${TOOLS_PATH}/genomebrowser/annotations ${BACKUPDIR_PATH}/
  mkdir -p ${TOOLS_PATH}/genomebrowser/annotations

  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Arabidopsis_lyrata.v.1.0.17.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Arabidopsis_thaliana.TAIR10.17.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Bos_taurus.UMD3.1.70.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Canis_familiaris.BROADD2.67.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Canis_familiaris.CanFam3.1.70.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Gallus_gallus.Gallus_gallus-4.0.pre.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Gasterosteus_aculeatus.BROADS1.70.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Homo_sapiens.GRCh37.70.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Homo_sapiens.NCBI36.54.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Human-MT.NC_012920.1.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Mus_musculus.GRCm38.70.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Mus_musculus.NCBIM37.67.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/N916Ysi.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Ovis_aries_v3.1.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/R1-RT.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Rattus_norvegicus.RGSC3.4.69.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Sus_scrofa.Sscrofa10.2.70.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/Vitis_vinifera.IGGP_12x.17.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/

  wget -O ${TOOLS_PATH}/genomebrowser/annotations/contents2.txt http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.0/contents2.txt

  echo "** Installing R library qvalue"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.12.1-vmbin/library/qvalue-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.12.1/lib64/R/library/

  echo "** Installing Rattus_norvegicus.Rnor_5.0.70"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Installing TagCleaner"
  curl -sL http://downloads.sourceforge.net/project/tagcleaner/standalone/tagcleaner-standalone-0.12.tar.gz | tar xz -C ${TOOLS_PATH}/
  ln -s tagcleaner-standalone-0.12 ${TOOLS_PATH}/tagcleaner

  echo "** Installing EMBOSS 6.5.7"
  #curl -s ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.5.7.tar.gz | tar -xz -C ${TMPDIR_PATH}/
  #cd ${TMPDIR_PATH}/EMBOSS-6.5.7
  #./configure --prefix=${TOOLS_PATH}/EMBOSS-6.5.7
  #make
  #make install
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/EMBOSS-6.5.7-vmbin.tar.gz | tar xz -C ${TOOLS_PATH}/
  ln -s EMBOSS-6.5.7 ${TOOLS_PATH}/emboss
  
  echo "** Installing fseq"
  curl -s http://fureylab.med.unc.edu/fseq/fseq_1.84.tgz | tar -xz -C ${TMPDIR_PATH}/ 
  mv -b ${TMPDIR_PATH}/fseq ${TOOLS_PATH}/fseq-1.84
  ln -s fseq-1.84 ${TOOLS_PATH}/fseq
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/misc/read_extend_bed.pm -o ${TOOLS_PATH}/fseq/bin/read_extend_bed.pm
  chmod 775 ${TOOLS_PATH}/fseq/bin/read_extend_bed.pm

  echo "** Installing R library zinba and dependencies"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1-vmbin/library/zinba.tar.gz | tar xz -C ${TOOLS_PATH}/R-2.15.1/lib64/R/library/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1-vmbin/library/R.utils.tar.gz | tar xz -C ${TOOLS_PATH}/R-2.15.1/lib64/R/library/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1-vmbin/library/doMC.tar.gz | tar xz -C ${TOOLS_PATH}/R-2.15.1/lib64/R/library/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1-vmbin/library/foreach.tar.gz | tar xz -C ${TOOLS_PATH}/R-2.15.1/lib64/R/library/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1-vmbin/library/multicore.tar.gz | tar xz -C ${TOOLS_PATH}/R-2.15.1/lib64/R/library/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1-vmbin/library/quantreg.tar.gz | tar xz -C ${TOOLS_PATH}/R-2.15.1/lib64/R/library/

  echo "** Installing rn5, mouse, human mirna indexes for bowtie, bowtie2"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_miRBase19_Rattus_norvegicus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_miRBase19_Rattus_norvegicus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_miRBase19_Homo_sapiens.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_miRBase19_Homo_sapiens.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_miRBase19_Mus_musculus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_miRBase19_Mus_musculus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
  mv -b ${TOOLS_PATH}/bowtie/indexes/mmu_miRB17mature.* ${BACKUPDIR_PATH}/

fi


# 2.5.1
compare_to_current_and_latest "2.5.1"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Installing rn5, mouse, human mirna indexes for bowtie, bowtie2"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_miRBase19_Rattus_norvegicus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_miRBase19_Rattus_norvegicus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_miRBase19_Homo_sapiens.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_miRBase19_Homo_sapiens.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_miRBase19_Mus_musculus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_miRBase19_Mus_musculus.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/

  # remove wrong target location for mirbase bowtie indexes in 2.5.0 first version
  rm -rf ${TOOLS_PATH}/indexes/

fi

# 2.5.2
compare_to_current_and_latest "2.5.2"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Installing genome browser fasta files"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Drosophila_melanogaster.BDGP5.70.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Updating genome browser annotations"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.2/Drosophila_melanogaster.BDGP5.70.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/

  wget -O ${TOOLS_PATH}/genomebrowser/annotations/contents2.txt http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.5.2/contents2.txt

  echo "** Updating HTSeq"
  ln -s /usr/local/bin/htseq-count_chr ${TOOLS_PATH}/htseq/htseq-count_chr

  sudo pip install HTSeq==0.5.4p3
  sudo wget -O /usr/local/bin/htseq-count_chr http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/htseq/htseq-count_chr 
  sudo chmod 755 /usr/local/bin/htseq-count_chr
  sudo wget -O /usr/local/lib/python2.7/dist-packages/HTSeq/scripts/count_chr.py http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/htseq/count_chr.py

  echo "** Installing announcements"
  sudo mkdir /opt/chipster-announcements
  sudo chown chipster /opt/chipster-announcements
  wget -O /opt/chipster-announcements/get-chipster-announcements.sh http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/virtual_machines/updates/misc/get-chipster-announcements.sh  
  sudo chmod 755 /opt/chipster-announcements/get-chipster-announcements.sh
  sudo wget -O /etc/update-motd.d/32-announcements http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/virtual_machines/updates/misc/32-announcements
  sudo chmod 755 /etc/update-motd.d/32-announcements      
  sudo wget -O /etc/cron.d/chipster-announcements http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/virtual_machines/updates/misc/chipster-announcements.crontab
  sudo chmod 644 /etc/cron.d/chipster-announcements

fi


# 2.6.0
compare_to_current_and_latest "2.6.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Installing png library for R-2.15.1_bioc-2.11"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1_bioc-2.11/library/png-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.15.1_bioc-2.11/lib64/R/library/
  # install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="png")

  echo "** Updating R-2.15.1"
  mv -b ${TOOLS_PATH}/R-2.15.1 ${BACKUPDIR_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1-vmbin/R-2.15.1-2013-05-26.tar.gz | tar -xz -C ${TOOLS_PATH}/

  #biocLite("crlmm")
  #biocLite("oligo")
  #biocLite("pd.mapping50k.hind240")
  #biocLite("pd.mapping50k.xba240")
  #biocLite("pd.mapping250k.nsp")
  #biocLite("pd.mapping250k.sty")
  #biocLite("pd.genomewidesnp.5")
  #biocLite("pd.genomewidesnp.6")
  #biocLite("genomewidesnp6Crlmm")
  #biocLite("genomewidesnp5Crlmm")
  #biocLite("human1mduov3bCrlmm") #Illumina
  #biocLite("human1mv1cCrlmm") #Illumina
  #biocLite("human370quadv3cCrlmm") #Illumina
  #biocLite("human370v1cCrlmm") #Illumina
  #biocLite("human550v3bCrlmm") #Illumina
  #biocLite("human610quadv1bCrlmm") #Illumina
  #biocLite("human650v3aCrlmm") #Illumina
  #biocLite("human660quadv1aCrlmm") #Illumina
  #biocLite("humanomni1quadv1bCrlmm") #Illumina
  #biocLite("humanomniexpress12v1bCrlmm") #Illumina
  

  echo "** Installing affy_20 for R-2.12"
  cd ${TMPDIR_PATH}/
  
  # TAIR annotations for Arabidosis 1.0 and 1.1
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/tairg.download/aragene11stattairgcdf_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/tairg.download/aragene11stattairgprobe_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/tairg.download/aragene10stattairgcdf_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/tairg.download/aragene10stattairgprobe_17.0.0.tar.gz'

  # Entrez annotations for human
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/hugene20sthsentrezgcdf_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/hugene20sthsentrezgprobe_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/hugene21sthsentrezgcdf_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/hugene21sthsentrezgprobe_17.0.0.tar.gz'

  # Entrez annotations for Arabidosis 1.0 and 1.1
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/aragene10statentrezgcdf_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/aragene10statentrezgprobe_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/aragene11statentrezgprobe_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/aragene11statentrezgcdf_17.0.0.tar.gz'

  # Entrez annotations for ragene_20 and ragene_21
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/ragene20strnentrezgcdf_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/ragene20strnentrezgprobe_17.0.0.tar.gz'
  ##wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/ragene20strnentrezg.db_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/ragene21strnentrezgprobe_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/ragene21strnentrezgcdf_17.0.0.tar.gz'

  # Entrez annotations for mouse_20 and mouse_21
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/mogene20stmmentrezgcdf_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/mogene20stmmentrezgprobe_17.0.0.tar.gz'
  ##wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/mogene20stmmentrezg.db_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/mogene21stmmentrezgcdf_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/mogene21stmmentrezgprobe_17.0.0.tar.gz'
  ##wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/mogene21stmmentrezg.db_17.0.0.tar.gz'

  # Entrez annotations for Danio rerio 1.0 and 1.1
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/zebgene10stdrentrezgcdf_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/zebgene10stdrentrezgprobe_17.0.0.tar.gz'
  ##wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/zebgene10stdrentrezg.db_17.0.0.tar.gz'

  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/zebgene11stdrentrezgprobe_17.0.0.tar.gz'
  ##wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/zebgene11stdrentrezg.db_17.0.0.tar.gz'
  #wget 'http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/17.0.0/entrezg.download/zebgene11stdrentrezgcdf_17.0.0.tar.gz'
  
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene11stattairgcdf_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene11stattairgprobe_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene10stattairgcdf_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene10stattairgprobe_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/hugene20sthsentrezgcdf_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/hugene20sthsentrezgprobe_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/hugene21sthsentrezgcdf_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/hugene21sthsentrezgprobe_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene10statentrezgcdf_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene10statentrezgprobe_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene11statentrezgprobe_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene11statentrezgcdf_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/ragene20strnentrezgcdf_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/ragene20strnentrezgprobe_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/ragene21strnentrezgprobe_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/ragene21strnentrezgcdf_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/mogene20stmmentrezgcdf_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/mogene20stmmentrezgprobe_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/mogene21stmmentrezgcdf_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/mogene21stmmentrezgprobe_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/zebgene10stdrentrezgcdf_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/zebgene10stdrentrezgprobe_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/zebgene11stdrentrezgprobe_17.0.0.tar.gz'
  wget 'http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/zebgene11stdrentrezgcdf_17.0.0.tar.gz'

  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ aragene10statentrezgcdf_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ aragene10statentrezgprobe_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ aragene10stattairgcdf_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ aragene10stattairgprobe_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ aragene11statentrezgcdf_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ aragene11statentrezgprobe_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ aragene11stattairgcdf_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ aragene11stattairgprobe_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ hugene20sthsentrezgcdf_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ hugene20sthsentrezgprobe_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ hugene21sthsentrezgcdf_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ hugene21sthsentrezgprobe_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ mogene20stmmentrezgcdf_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ mogene20stmmentrezgprobe_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ mogene21stmmentrezgcdf_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ mogene21stmmentrezgprobe_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ ragene20strnentrezgcdf_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ ragene20strnentrezgprobe_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ ragene21strnentrezgcdf_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ ragene21strnentrezgprobe_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ zebgene10stdrentrezgcdf_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ zebgene10stdrentrezgprobe_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ zebgene11stdrentrezgcdf_17.0.0.tar.gz
  ${TOOLS_PATH}/R-2.12.1/bin/R CMD INSTALL -l ${TOOLS_PATH}/R-2.12.1/lib64/R/library/ zebgene11stdrentrezgprobe_17.0.0.tar.gz
  
  
  echo "** Add fasta links to bowtie2 indexes"
  rm -f ${TOOLS_PATH}/bowtie2/indexes/Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa
  ln -s ../../genomes/fasta/nochr/Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa ${TOOLS_PATH}/bowtie2/indexes
  rm -f ${TOOLS_PATH}/bowtie/indexes/Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa
  ln -s ../../genomes/fasta/nochr/Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa ${TOOLS_PATH}/bowtie/indexes

  echo "** Updating Sus_scrofa bowtie2 index"
  rm ${TOOLS_PATH}/bowtie2/indexes/Sus_scrofa.Sscrofa10.2.69.dna.toplevel.*
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Sus_scrofa.Sscrofa10.2.69.tar.gz | tar -xz -C ${TOOLS_PATH}/bowtie2/
        
        
fi



# 2.6.1
compare_to_current_and_latest "2.6.1"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Fixing HTSeq"
  rm -f ${TOOLS_PATH}/htseq/htseq-count_chr
  rm -f ${TOOLS_PATH}/htseq/htseq_count_chr
  ln -s /usr/local/bin/htseq-count_chr ${TOOLS_PATH}/htseq/htseq-count_chr
  
  echo "** Updating genomes"
  rm -f ${TOOLS_PATH}/bowtie/indexes/canFam2*
  rm -f ${TOOLS_PATH}/bowtie2/indexes/canFam3*

  rm -f ${TOOLS_PATH}/bowtie/indexes/e_coli*
  rm -f ${TOOLS_PATH}/bowtie2/indexes/e_coli*
  rm -f ${TOOLS_PATH}/genomes/fasta/e_coli.fa

  rm -f ${TOOLS_PATH}/bowtie/indexes/miRBase18_mmu_matureT.fa*
  rm -f ${TOOLS_PATH}/bwa_indexes/mmu_miRB17mature*

  rm -f ${TOOLS_PATH}/bowtie/indexes/Gasterosteus_aculeatus.BROADS1.67* 
  
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_athaliana.TAIR10.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_athaliana.TAIR10.tar.gz | tar xz -C ${TOOLS_PATH}/
  
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_Canis_familiaris.CanFam3.1.71.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Canis_familiaris.CanFam3.1.71.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_Canis_familiaris.CanFam3.1.71.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Canis_familiaris.CanFam3.1.71.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Canis_familiaris.CanFam3.1.71.tar.gz | tar xz -C ${TOOLS_PATH}/
  
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_Escherichia_coli_n1.GCA_000303635.1.18.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Escherichia_coli_n1.GCA_000303635.1.18.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_Escherichia_coli_n1.GCA_000303635.1.18.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Escherichia_coli_n1.GCA_000303635.1.18.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Escherichia_coli_n1.GCA_000303635.1.18.tar.gz | tar xz -C ${TOOLS_PATH}/
  
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_Gasterosteus_aculeatus.BROADS1.71.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie2_index_Gasterosteus_aculeatus.BROADS1.71.tar.gz | tar xz -C ${TOOLS_PATH}/bowtie2/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_Gasterosteus_aculeatus.BROADS1.71.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/fasta_nochr_Gasterosteus_aculeatus.BROADS1.71.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Gasterosteus_aculeatus.BROADS1.71.tar.gz | tar xz -C ${TOOLS_PATH}/
  
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Rattus_norvegicus.Rnor_5.0.70.tar.gz | tar xz -C ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Sus_scrofa.Sscrofa10.2.69.tar.gz | tar xz -C ${TOOLS_PATH}/
  
  echo "** Updating cufflinks2"
  curl -s http://cufflinks.cbcb.umd.edu/downloads/cufflinks-2.1.1.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/
  rm -f ${TOOLS_PATH}/cufflinks2
  ln -s cufflinks-2.1.1.Linux_x86_64 ${TOOLS_PATH}/cufflinks2
  rm -rf ${TOOLS_PATH}/cufflinks-2.0.2.Linux_x86_64
              
fi

# 2.7.0
compare_to_current_and_latest "2.7.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Installing sva R library for R-2.15_bioc-2.11"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1_bioc-2.11/library/sva-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.15.1_bioc-2.11/lib64/R/library/
  
  echo "** Installing mothur"
  cd ${TMPDIR_PATH}/
  wget -nv http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/mothur/Mothur-1.28.cen_64.noReadLine.zip
  unzip -q Mothur-1.28.cen_64.noReadLine.zip
  mv mothur ${TOOLS_PATH}/mothur-1.28
  ln -s mothur-1.28 ${TOOLS_PATH}/mothur
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/mothur/mothur-data.tar.gz | tar -xz -C ${TOOLS_PATH}/
                                                                                                                                                                                                                                                                                                        
fi

# 2.7.2
compare_to_current_and_latest "2.7.2"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  echo "** Updating Rattus_norvegicus.RGSC3.4.68.gtf"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/genomes/gtf_Rattus_norvegicus.RGSC3.4.68.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Updating R-2.15 with Bioconductor 2.11"
  # biocLite("WGCNA")
  mv -b ${TOOLS_PATH}/R-2.15.1_bioc-2.11 ${BACKUPDIR_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1_bioc-2.11/R-2.15.1_bioc-2.11-vmbin_v4.tar.gz | tar -xz -C ${TOOLS_PATH}/  

  echo "** Installing Mfuzz and dependencies for R-2.12"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.12.1-vmbin/library/widgetTools-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.12.1/lib64/R/library/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.12.1-vmbin/library/tkWidgets-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.12.1/lib64/R/library/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.12.1-vmbin/library/Mfuzz-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.12.1/lib64/R/library/

fi



# 2.8.0
compare_to_current_and_latest "2.8.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  # Update Bowtie 2
  echo "** Updating Bowtie 2"
  cd ${TMPDIR_PATH}/
  wget -nv -O bowtie2-2.1.0-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip/download
  unzip -q bowtie2-2.1.0-linux-x86_64.zip
  mv bowtie2-2.1.0 ${TOOLS_PATH}
  mv ${TOOLS_PATH}/bowtie2/indexes ${TOOLS_PATH}/bowtie2-2.1.0/ 
  mv ${TOOLS_PATH}/bowtie2-2.0.0-beta7 ${BACKUPDIR_PATH}/
  rm ${TOOLS_PATH}/bowtie2
  ln -s bowtie2-2.1.0 ${TOOLS_PATH}/bowtie2
  rm ${TMPDIR_PATH}/bowtie2-2.1.0-linux-x86_64.zip
 
  # Add RmiR.Hs.miRNA to R-2.15.1_bioc-2.11
  echo "** Installing RmiR.Hs.miRNA to R-2.15.1_bioc-2.11"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1_bioc-2.11/library/RmiR.Hs.miRNA-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.15.1_bioc-2.11/lib64/R/library/

fi


compare_to_current_and_latest "2.8.1"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 

  # Update Tophat 2
  echo "** Updating Tophat 2"
  curl -s http://tophat.cbcb.umd.edu/downloads/tophat-2.0.9.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/
  mv ${TOOLS_PATH}/tophat-2.0.4.Linux_x86_64 ${BACKUPDIR_PATH}/
  rm ${TOOLS_PATH}/tophat2
  ln -s tophat-2.0.9.Linux_x86_64 ${TOOLS_PATH}/tophat2

  # ** Install tabix
  echo "** Installing tabix"
  # curl -sL http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2/download | tar xj -C ${TOOLS_PATH}/
  # cd ${TOOLS_PATH}/tabix-0.2.6
  # make
  # ln -s tabix-0.2.6 ${TOOLS_PATH}/tabix 
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/tabix-0.2.6-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s tabix-0.2.6 ${TOOLS_PATH}/tabix 

fi

# 2.9.0
compare_to_current_and_latest "2.9.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then
 
  echo "** Installing genome bundle tool dependencies"
  sudo apt-get -y install python3-yaml
  touch installed.yaml

  # bundle tool will replace these with newer versions
  mv -b ${TOOLS_PATH}/genomebrowser/annotations/Drosophila_melanogaster.BDGP5.70* ${BACKUPDIR_PATH}/
  mv -b ${TOOLS_PATH}/genomes/fasta/nochr/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa ${BACKUPDIR_PATH}/

  echo "** Removing FREEC_Linux64"
  rm -rf ${TOOLS_PATH}/FREEC_Linux64

  echo "** Installing R-3.0"
  curl -L http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.0.2-vmbin/R-3.0.2-2013-11-22.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Updating bedtools"
  mv -b ${TOOLS_PATH}/BEDTools-Version-2.12.0 ${BACKUPDIR_PATH}/
  #cd ${TMPDIR_PATH}/
  #curl -s http://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz | tar -xz
  #cd bedtools-2.17.0
  #make clean
  #make all
  #cd ../
  #mv bedtools-2.17.0 ${TOOLS_PATH}/
  curl -L http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bedtools-2.17.0-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  rm ${TOOLS_PATH}/bedtools
  ln -s bedtools-2.17.0 ${TOOLS_PATH}/bedtools

  echo "** Updating samtools"
  mv -b ${TOOLS_PATH}/samtools-0.1.18 ${BACKUPDIR_PATH}/
  #cd ${TMPDIR_PATH}/
  #curl -sL http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download | tar -xj
  #cd samtools-0.1.19/
  #make
  #cd ../
  #mv samtools-0.1.19/ ${TOOLS_PATH}
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/samtools-0.1.19-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  rm ${TOOLS_PATH}/samtools
  ln -s samtools-0.1.19 ${TOOLS_PATH}/samtools
  
  echo "** Updating vcftools"
  mv -b ${TOOLS_PATH}/vcftools_0.1.9 ${BACKUPDIR_PATH}/
  #cd ${TMPDIR_PATH}/
  #curl -sL http://sourceforge.net/projects/vcftools/files/vcftools_0.1.11.tar.gz/download| tar -xz
  #cd vcftools_0.1.11/
  #make
  #cd ../
  #mv vcftools_0.1.11/ ${TOOLS_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/vcftools_0.1.11-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  rm ${TOOLS_PATH}/vcftools
  ln -s vcftools_0.1.11 ${TOOLS_PATH}/vcftools
            
  echo "** Installing ConsensuPathDB tool dependency"
  sudo apt-get -y install python-zsi
fi

# 2.11.0
compare_to_current_and_latest "2.11.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then
  
  echo "** Installing Picard tools"
  cd ${TOOLS_PATH}
  wget -nv -O picard-tools-1.105.zip http://sourceforge.net/projects/picard/files/picard-tools/1.105/picard-tools-1.105.zip/download
  unzip -q picard-tools-1.105.zip
  ln -s picard-tools-1.105 picard-tools
  # remove this optional package because it's in the root of the tools
  rm snappy-java-1.0.3-rc3.jar
  rm picard-tools-1.105.zip

  echo "** Installing RSeQC"
  cd ${TOOLS_PATH}
  curl -L http://sourceforge.net/projects/rseqc/files/RSeQC-2.3.7.tar.gz/download | tar -xz
  ln -s RSeQC-2.3.7 ${TOOLS_PATH}/RSeQC
  cd RSeQC
  sudo python setup.py install

  echo "** Installing EMBOSS"
  # dependency for png images
  #sudo apt-get -y install libgd2-noxpm-dev
  # update EMBOSS
  mv -b ${TOOLS_PATH}/EMBOSS-6.5.7 ${BACKUPDIR_PATH}/
  curl http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/EMBOSS-6.5.7-with-extras-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Installing primer3"
  curl http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/primer3-vmbin.tar.gz| tar -xz -C ${TOOLS_PATH}/

  echo "** Installing meme"
  curl http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/meme_4.2.0-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s meme_4.2.0 ${TOOLS_PATH}/meme

  echo "** Installing human tophat index"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/misc/hg19.ti.tar.gz | tar -xzv -C ${TOOLS_PATH}/bowtie2/indexes/

  echo "** Installing GRCh37_74 ensembl transcripts"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/misc/GRCh37_74_ensembl_transcripts.tar.gz | tar -xzv -C ${TOOLS_PATH}/bowtie2/indexes/

  echo "** Installing dimont"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/misc/dimont.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Installing Trimmomatic"
  cd ${TMPDIR_PATH}/
  wget -nv http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
  unzip Trimmomatic-0.32.zip
  mv Trimmomatic-0.32 ${TOOLS_PATH}/
  ln -s Trimmomatic-0.32 ${TOOLS_PATH}/trimmomatic

  echo "** Installing express"
  curl http://bio.math.berkeley.edu/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz | tar -xz -C ${TOOLS_PATH}/
  ln -s express-1.5.1-linux_x86_64 ${TOOLS_PATH}/express
 
  echo "** Installing blast"       
  curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.29+-x64-linux.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s ncbi-blast-2.2.29+ ${TOOLS_PATH}/blast                                  

  echo "** Updating prinseq"
  cd ${TMPDIR_PATH}/
  curl -L http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz/download | tar -xz
  chmod a+x prinseq-lite-0.20.4/*.pl
  mv -b ${TOOLS_PATH}/prinseq-lite-0.19.3 ${BACKUPDIR_PATH}/
  mv prinseq-lite-0.20.4 ${TOOLS_PATH}/
  rm ${TOOLS_PATH}/prinseq
  ln -s prinseq-lite-0.20.4 ${TOOLS_PATH}/prinseq
  
  echo "** Installing mafft"
  #cd ${TMPDIR_PATH}/
  #curl http://mafft.cbrc.jp/alignment/software/mafft-7.130-without-extensions-src.tgz | tar -xz
  #cd mafft-7.130-without-extensions/core
  #sed -i 's/PREFIX = \/usr\/local/PREFIX = \/opt\/chipster\/tools\/mafft-7.130-without-extensions/g' Makefile
  #make clean
  #make
  #make install
  curl http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/misc/mafft-7.130-without-extensions-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/
  ln -s mafft-7.130-without-extensions ${TOOLS_PATH}/mafft

  echo "** Updating R-3.0.2"
  rm -rf ${TOOLS_PATH}/R-3.0.2
  curl -L http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.0.2-vmbin/R-3.0.2-2014-03-03.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Updating DEXSeq scripts"
  cd ${TMPDIR_PATH}/
  curl -sL http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/DEXSeq_1.8.0.tar.gz | tar -xz
  cp DEXSeq/inst/python_scripts/dexseq_count.py ${TOOLS_PATH}/dexseq-exoncounts  
  cp DEXSeq/inst/python_scripts/dexseq_prepare_annotation.py ${TOOLS_PATH}/dexseq-exoncounts

  # Update Tophat 2
  echo "** Updating Tophat 2"
  curl -s http://tophat.cbcb.umd.edu/downloads/tophat-2.0.10.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/
  mv ${TOOLS_PATH}/tophat-2.0.9.Linux_x86_64 ${BACKUPDIR_PATH}/
  rm ${TOOLS_PATH}/tophat2
  ln -s tophat-2.0.10.Linux_x86_64 ${TOOLS_PATH}/tophat2
fi

# 2.11.1
compare_to_current_and_latest "2.11.1"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then

  echo "** Installing R libs"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.0.2-vmbin/library/png-0.1-7-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-3.0.2/lib64/R/library/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.0.2-vmbin/library/QDNAseq-0.99.4-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-3.0.2/lib64/R/library/
  rm -rf ${TOOLS_PATH}/R-3.0.2/lib64/R/library/WECCA
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-3.0.2-vmbin/library/WECCA-0.40-vmbin.tar.gz | tar -xz -C ${TOOLS_PATH}/R-3.0.2/lib64/R/library/

  echo "** Installing QDNAseq files"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/misc/QDNAseq.hg19.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Removing obsolete files"
  rm -rf ${TOOLS_PATH}/MPScall
  rm -rf ${TOOLS_PATH}/CanGEM
  rm -rf ${TOOLS_PATH}/DGV

fi

# 2.12.0
compare_to_current_and_latest "2.12.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then

  echo "** Updating genome browser annotation index"
  # remove old droso
  wget -O ${TOOLS_PATH}/genomebrowser/annotations/contents2.txt http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/2.11.2/contents2.txt

fi

# 2.12.1
compare_to_current_and_latest "2.12.1"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then
  echo ""
fi



#####################################
# VERSION SPECIFIC ENTRIES END HERE #
#####################################

# Update Chipster itself (incl. tool scripts), unless already at latest
compare_to_current "$LATEST_VERSION"
if [ $CURRENT_COMPARED -lt 0 ] ; then 

  echo "** Updating Chipster installation to $LATEST_VERSION"
  cd ${CHIP_PATH}/
  
  # Get install package (override, if exists)
  rm -f chipster-$LATEST_VERSION.tar.gz
    wget http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/versions/$LATEST_VERSION/chipster-$LATEST_VERSION.tar.gz

  # Unpack libs
    echo "** Updating Chipster libs: shared/libs"
    mv shared ${BACKUPDIR_PATH}/
    tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/shared

  # Unpack webstat web-root including client jar
    echo "** Updating Chipster web: webstart/web-root"
    cp webstart/web-root/chipster.jnlp ${BACKUPDIR_PATH}/
    cp webstart/web-root/chipster-config.xml ${BACKUPDIR_PATH}/
    mv webstart/web-root ${BACKUPDIR_PATH}/
    tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/webstart/web-root
    cp ${BACKUPDIR_PATH}/chipster.jnlp webstart/web-root/ 
    cp ${BACKUPDIR_PATH}/chipster-config.xml webstart/web-root/ 

  # Copy away tool scripts in case there were important local changes
    cp -r comp/modules ${BACKUPDIR_PATH}/

  # Unpack bundle tool
    echo "** Updating Chipster genome bundle tool"
    tar -C .. --overwrite -xzf chipster-$LATEST_VERSION.tar.gz chipster/bundle.py

  # Unpack tool scripts
    echo "** Updating Chipster tool scripts: comp/modules"
    tar -C .. --overwrite -xzf chipster-$LATEST_VERSION.tar.gz chipster/comp/modules

  # Update runtimes.xml
  echo "** Updating Chipster runtimes: comp/conf/runtimes.xml"
  cp -r comp/conf/runtimes.xml ${BACKUPDIR_PATH}/
  tar -C .. --overwrite -xzf chipster-$LATEST_VERSION.tar.gz chipster/comp/conf/runtimes.xml


  # Update webapps
  rm -rf webstart/webapps
  tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/webstart/webapps/tool-editor.war
  rm -rf manager/webapps
  tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/manager/webapps/admin-web.war

  # Update R libs
  #echo "** Updating R libraries"
  #${TOOLS_PATH}/R-${R_VERSION}/bin/Rscript --vanilla ${CHIP_PATH}/comp/modules/admin/R/install-libs.R

  # Clean up
  rm chipster-$LATEST_VERSION.tar.gz
fi

# Remove temp dir
rm -rf ${TMPDIR_PATH}/

# Bundle
function update_bundles()
{
  echo "** Updating genome bundles"
  wget -q http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bundle/bundles.yaml -O bundles.yaml
  python3 bundle.py update installed -q
}

function install_bundle()
{
  echo "** Installing $1 genome"
  python3 bundle.py install $1.bowtie
  python3 bundle.py install $1.bowtie2
  python3 bundle.py install $1.bwa
  python3 bundle.py install $1.gb
  python3 bundle.py install $1
}

# Version specific bundle tool commands
compare_to_current_and_latest "2.9.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then 
  update_bundles
  install_bundle "Drosophila_melanogaster.BDGP5"
fi

# 2.12.0
compare_to_current_and_latest "2.12.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then
  update_bundles
  install_bundle "Schizosaccharomyces_pombe.ASM294v2"
fi

# Check backup dir
SIZE=`du -hs ${BACKUPDIR_PATH} | cut -f1`
echo ""
echo "Total of $SIZE old data has been backed up to ${BACKUPDIR_PATH}"
echo "It is recommended to inspect the directory and then to remove it"
   
# We are done
echo "Update completed successfully"
echo "Remember to start the Chipster service: 'service chipster start'"
echo $END_MESSAGE
