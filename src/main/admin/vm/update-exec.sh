#!/bin/bash

# This script updates to latest version. Updates between minor versions should be smooth and
# automatic, where as updates between major versions can require some manual steps afterwards
# if some specific local customisations were in place.
# This update mechanism has been available since 2.0.2.

# Latest version, matching tar-packages must be available 
LATEST_VERSION=2.1.0

# Exit immediately if some command fails
set -e

# Detect current version
CURRENT_VERSION=`ls -1 shared/lib | grep ^chipster-[0-9\\.]*.jar | gawk 'match($0, "chipster-([0-9\\\\.]*).jar", g) {print g[1]}'`
CURRENT_MAIN_VERSION=`echo $CURRENT_VERSION | gawk 'match($0, "([0-9]*).[0-9]*.[0-9]*", g) {print g[1]}'`
CURRENT_MAJOR_VERSION=`echo $CURRENT_VERSION | gawk 'match($0, "[0-9]*.([0-9]*).[0-9]*", g) {print g[1]}'`
CURRENT_MINOR_VERSION=`echo $CURRENT_VERSION | gawk 'match($0, "[0-9]*.[0-9]*.([0-9]*)", g) {print g[1]}'`
LATEST_MAIN_VERSION=`echo $LATEST_VERSION | gawk 'match($0, "([0-9]*).[0-9]*.[0-9]*", g) {print g[1]}'`
LATEST_MAJOR_VERSION=`echo $LATEST_VERSION | gawk 'match($0, "[0-9]*.([0-9]*).[0-9]*", g) {print g[1]}'`
LATEST_MINOR_VERSION=`echo $LATEST_VERSION | gawk 'match($0, "[0-9]*.[0-9]*.([0-9]*)", g) {print g[1]}'`

# Check if versions match
echo Detected version $CURRENT_VERSION
if [ $CURRENT_MAIN_VERSION -ge $LATEST_MAIN_VERSION -a $CURRENT_MAJOR_VERSION -ge $LATEST_MAJOR_VERSION -a $CURRENT_MINOR_VERSION -ge $LATEST_MINOR_VERSION ]; then
	echo "Already at latest version, nothing needs to be updated"
	exit
fi

# Confirm update
echo "Will update to version $LATEST_VERSION"
echo "Update will start next. It can take several hours, depending on your network connection"
echo "Do you wish to proceed?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) echo "Update started"; break;;
        No ) echo "Update aborted"; exit;;
    esac
done

# Start update
INST_PATH=/opt
CHIP_PATH=${INST_PATH}/chipster
TOOLS_PATH=${CHIP_PATH}/tools
TMPDIR_PATH=/tmp/chipster-install-temp

# Create temp dir
rm -rf ${TMPDIR_PATH}/
mkdir ${TMPDIR_PATH}/

# Create backup dir
TIMESTAMP=$(date +"%Y%m%d%H%M")
BACKUPDIR_PATH=/tmp/$TIMESTAMP
mkdir ${BACKUPDIR_PATH}/

#######################################
# VERSION SPECIFIC ENTRIES START HERE #
# (ADD NEW ENTRIES TO THE END)        #
#######################################

# 2.0.3
if [ $CURRENT_MAIN_VERSION -lt 2 -o  $CURRENT_MAJOR_VERSION -lt 0 -o $CURRENT_MINOR_VERSION -lt 3 ] ; then
	echo "Installing mm10 bowtie indexes"
	curl -L http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_mm10.tar.gz  | tar -xz -C ${TOOLS_PATH}/bowtie/
fi

# 2.1.0 
if [ $CURRENT_MAIN_VERSION -lt 2 -o  $CURRENT_MAJOR_VERSION -lt 1 -o $CURRENT_MINOR_VERSION -lt 0 ] ; then

    echo "Updating prinseq"
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


# 2.1.1 
if [ $CURRENT_MAIN_VERSION -lt 2 -o  $CURRENT_MAJOR_VERSION -lt 1 -o $CURRENT_MINOR_VERSION -lt 1 ] ; then

		echo "Updating mm10 bowtie index files
  	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_mm10.tar.gz  | tar -xz -C ${TOOLS_PATH}/bowtie/

		echo "Updating samtools
		cd ${TMPDIR_PATH}/
  	curl -sL http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2/download | tar -xj
  	cd samtools-0.1.18/
  	make
  	cd ../
  	mv samtools-0.1.18/ ${TOOLS_PATH}
  	mv ${TOOLS_PATH}/samtools-0.1.13/ ${BACKUPDIR_PATH}/
  	rm ${TOOLS_PATH}/samtools
  	ln -s samtools-0.1.18 ${TOOLS_PATH}/samtools

		echo "Installing bowtie2"
  	cd ${TMPDIR_PATH}/
  	wget -nv -O bowtie2-2.0.0-beta7-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.0.0-beta7/bowtie2-2.0.0-beta7-linux-x86_64.zip/download
  	unzip -q bowtie2-2.0.0-beta7-linux-x86_64.zip
  	mv bowtie2-2.0.0-beta7 ${TOOLS_PATH}
  	ln -s bowtie2-2.0.0-beta7 ${TOOLS_PATH}/bowtie2
  	rm bowtie2-2.0.0-beta7-linux-x86_64.zip

		echo "Installing bowtie2 indexes"
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
	
		echo "Installing tophat2
  	cd ${TMPDIR_PATH}/
  	curl -s http://tophat.cbcb.umd.edu/downloads/tophat-2.0.4.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/
  	ln -s tophat-2.0.4.Linux_x86_64 ${TOOLS_PATH}/tophat2

    echo "Installing DEXSeq"
    # DEXSeq
		cd ${TMPDIR_PATH}/
		curl -sL http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/DEXSeq_1.2.1.tar.gz | tar -xz
		mkdir ${TOOLS_PATH}/dexseq-exoncounts
		cp DEXSeq/inst/python_scripts/dexseq_count.py ${TOOLS_PATH}/dexseq-exoncounts  
		cp DEXSeq/inst/python_scripts/dexseq_prepare_annotation.py ${TOOLS_PATH}/dexseq-exoncounts
		rm -rf DEXSeq	

	echo "Install new cufflinks" 
	cd ${TMPDIR_PATH}/
	curl -s http://cufflinks.cbcb.umd.edu/downloads/cufflinks-2.0.2.Linux_x86_64.tar.gz | tar -xz -C ${TOOLS_PATH}/
                        
fi



#####################################
# VERSION SPECIFIC ENTRIES END HERE #
#####################################

# Update Chipster itself (incl. tool scripts), unless already at latest
if [ $CURRENT_MAIN_VERSION -lt $LATEST_MAIN_VERSION -o  $CURRENT_MAJOR_VERSION -lt $LATEST_MAJOR_VERSION -o $CURRENT_MINOR_VERSION -lt $LATEST_MINOR_VERSION ] ; then

	echo "Updating Chipster installation to $LATEST_VERSION"
	cd ${CHIP_PATH}/
	
	# Get install package (override, if exists)
	rm -f chipster-$LATEST_VERSION.tar.gz
    wget http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/versions/$LATEST_VERSION/chipster-$LATEST_VERSION.tar.gz

	# Move away old libs to avoid conflicts when lib names change
    mv shared ${BACKUPDIR_PATH}/
    mv webstart/web-root/lib ${BACKUPDIR_PATH}/

	# Unpack libs
    echo "Updating libs: shared/libs"
    tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/shared
    echo "Updating libs: webstart/web-root/lib"
    tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/webstart/web-root/lib

	# Copy away tool scripts in case there were important local changes
    cp -r comp/modules ${BACKUPDIR_PATH}/

	# Unpack tool scripts
    echo "Updating tool scripts: comp/modules"
    tar -C .. --overwrite -xzf chipster-$LATEST_VERSION.tar.gz chipster/comp/modules

	# Update manuals
	echo "Updating manuals"
	mv webstart/web-root/manual ${BACKUPDIR_PATH}/
	tar -C .. --overwrite -xzf chipster-$LATEST_VERSION.tar.gz chipster/webstart/web-root/manual

	# Update runtimes.xml
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
