#!/bin/bash

# This script updates to latest version. Updates between minor versions should be smooth and
# automatic, where as updates between major versions can require some manual steps afterwards
# if some specific local customisations were in place.
# This update mechanism has been available since 2.0.2.

# Latest version, matching tar-packages must be available 
LATEST_VERSION=2.2.2

# Exit immediately if some command fails
set -e

# Helper function
function compare_to_current()
{
  VERSION_ARR=( `echo $1 | tr "." "\n"` )
  CURRENT_VERSION_ARR=( `echo $CURRENT_VERSION | tr "." "\n"`)

  # Check main version  
  if [ ${VERSION_ARR[0]} -lt ${CURRENT_VERSION_ARR[0]} ] ; then
    CURRENT_COMPARED=1
    return
  fi
  if [ ${VERSION_ARR[0]} -gt ${CURRENT_VERSION_ARR[0]} ] ; then
    CURRENT_COMPARED=-1
    return
  fi

  # Check major version  
  if [ ${VERSION_ARR[1]} -lt ${CURRENT_VERSION_ARR[1]} ] ; then
    CURRENT_COMPARED=1
    return
  fi
  if [ ${VERSION_ARR[1]} -gt ${CURRENT_VERSION_ARR[1]} ] ; then
    CURRENT_COMPARED=-1
    return
  fi
  
  # Check minor version
  if [ ${VERSION_ARR[2]} -lt ${CURRENT_VERSION_ARR[2]} ] ; then
    CURRENT_COMPARED=1
    return
  fi
  if [ ${VERSION_ARR[2]} -gt ${CURRENT_VERSION_ARR[2]} ] ; then
    CURRENT_COMPARED=-1
    return
  fi
  
  CURRENT_COMPARED=0
}

# Detect current version
CURRENT_VERSION=`ls -1 shared/lib | grep ^chipster-[0-9\\.]*.jar | gawk 'match($0, "chipster-([0-9\\\\.]*).jar", g) {print g[1]}'`

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

# Confirm update
echo "Will update to version $LATEST_VERSION"
echo "Update will start next. It can take several hours, depending on your network connection"
echo "Do you wish to proceed?"
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
compare_to_current "2.0.3"
if [ $CURRENT_COMPARED -lt 0 ] ; then 

	echo "** Installing mm10 bowtie indexes"
	curl -L http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_mm10.tar.gz  | tar -xz -C ${TOOLS_PATH}/bowtie/
fi

# 2.1.0 
compare_to_current "2.1.0"
if [ $CURRENT_COMPARED -lt 0 ] ; then 

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
    mv ${TOOLS_PATH}/genomebrowser/annotations ${BACKUPDIR_PATH}/
    mkdir ${TOOLS_PATH}/genomebrowser/annotations # not typically needed, but the tar package is a bit stupid in this case
    curl -L http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/All_genomes_for_browser_v2.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  
    echo "Installing R-2.15"
    curl -L http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1.tar.gz | tar -xz -C ${TOOLS_PATH}/
        
fi


# 2.2.0 
compare_to_current "2.2.0"
if [ $CURRENT_COMPARED -lt 0 ] ; then 

	echo "** Updating samtools"
	cd ${TMPDIR_PATH}/
  mv ${TOOLS_PATH}/samtools-0.1.13/ ${BACKUPDIR_PATH}/
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
	mv ${TOOLS_PATH}/genomes/Homo_sapiens.GRCh37.62.chr.gtf ${BACKUPDIR_PATH}/
	mv ${TOOLS_PATH}/genomes/Mus_musculus.NCBIM37.62.chr.gtf ${BACKUPDIR_PATH}/
	mv ${TOOLS_PATH}/genomes/Rattus_norvegicus.RGSC3.4.62.chr.gtf ${BACKUPDIR_PATH}/
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
compare_to_current "2.2.1"
if [ $CURRENT_COMPARED -lt 0 ] ; then 
  echo "** Tools image is uptodate"                      
fi

# 2.2.2
compare_to_current "2.2.2"
if [ $CURRENT_COMPARED -lt 0 ] ; then 
  echo "** Removing obsolete link and tophat"                      
  rm -f ${TOOLS_PATH}/admin/ngs
  rm -rf ${TOOLS_PATH}/tophat-1.3.0.Linux_x86_64	
  
  echo "** Adding sheep to genome"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/genomebrowser_fasta_Ovis_aries.Oar_v3.1.dna.toplevel.tar.gz | tar -xz -C ${TOOLS_PATH}/genomebrowser/annotations/
  wget -O ${TOOLS_PATH}/genomebrowser/annotations/contents2.txt http://www.nic.funet.fi/pub/sci/molbio/chipster/annotations/compressed/contents2.txt

  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bwa_indexes/bwa_index_ovis_aries_texel.tar.gz | tar -xz -C ${TOOLS_PATH}/
  
  echo "** Installing R library maSigPro"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/maSigPro.tar.gz | tar -xz -C ${TOOLS_PATH}/R-2.12.1/lib64/R/library/
  
  echo "** Updating R-2.15"
  mv ${TOOLS_PATH}/R-2.15.1 ${BACKUPDIR_PATH}/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R/R-2.15.1.tar.gz | tar -xz -C ${TOOLS_PATH}/

  echo "** Installing GATK"
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/GenomeAnalysisTKLite-latest.tar.bz2 | tar -xj -C ${TOOLS_PATH}/
  ln -s GenomeAnalysisTKLite-2.1-11-gfb37f33 ${TOOLS_PATH}/gatk

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

	# Move away old libs to avoid conflicts when lib names change
    mv shared ${BACKUPDIR_PATH}/
    mv webstart/web-root/lib ${BACKUPDIR_PATH}/

	# Unpack libs
    echo "** Updating Chipster libs: shared/libs"
    tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/shared
    echo "** Updating Chipster libs: webstart/web-root/lib"
    tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/webstart/web-root/lib

	# Copy away tool scripts in case there were important local changes
    cp -r comp/modules ${BACKUPDIR_PATH}/

	# Unpack tool scripts
    echo "** Updating Chipster tool scripts: comp/modules"
    tar -C .. --overwrite -xzf chipster-$LATEST_VERSION.tar.gz chipster/comp/modules

	# Update manuals
	echo "** Updating Chipster manuals: webstart/web-root/manual"
	mv webstart/web-root/manual ${BACKUPDIR_PATH}/
	tar -C .. --overwrite -xzf chipster-$LATEST_VERSION.tar.gz chipster/webstart/web-root/manual

	# Update runtimes.xml
	echo "** Updating Chipster runtimes: comp/conf/runtimes.xml"
    cp -r comp/conf/runtimes.xml ${BACKUPDIR_PATH}/
	tar -C .. --overwrite -xzf chipster-$LATEST_VERSION.tar.gz chipster/comp/conf/runtimes.xml

	# Clean up
    rm chipster-$LATEST_VERSION.tar.gz
fi

# Remove temp dir
rm -rf ${TMPDIR_PATH}/

# Check backup dir
SIZE=`du -hs ${BACKUPDIR_PATH} | cut -f1`
echo "Total of $SIZE old data has been backed up to ${BACKUPDIR_PATH}"
echo "It is recommended to inspect the directory and then to remove it"
   
# We are done
echo "Update completed successfully"
