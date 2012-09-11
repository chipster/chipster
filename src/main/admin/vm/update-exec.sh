#!/bin/bash

# This script updates to latest minor version of the same major version (e.g. 2.0.1 -> 2.0.3)
LATEST_VERSION=2.0.4

# Detect current version
CURRENT_VERSION=`ls -1 shared/lib | grep ^chipster-[0-9\\.]*.jar | gawk 'match($0, "chipster-([0-9\\\\.]*).jar", g) {print g[1]}'`
CURRENT_MAJOR_VERSION=`echo $CURRENT_VERSION | gawk 'match($0, "([0-9]*.[0-9]*).[0-9]*", g) {print g[1]}'`
CURRENT_MINOR_VERSION=`echo $CURRENT_VERSION | gawk 'match($0, "[0-9]*.[0-9]*.([0-9]*)", g) {print g[1]}'`
LATEST_MAJOR_VERSION=`echo $LATEST_VERSION | gawk 'match($0, "([0-9]*.[0-9]*).[0-9]*", g) {print g[1]}'`
LATEST_MINOR_VERSION=`echo $LATEST_VERSION | gawk 'match($0, "[0-9]*.[0-9]*.([0-9]*)", g) {print g[1]}'`


# Check if versions match
echo Detected version $CURRENT_VERSION
if [ ! "$CURRENT_MAJOR_VERSION" = "$LATEST_MAJOR_VERSION" ]; then
	echo "Error! Will work only for major version $LATEST_MAJOR_VERSION!"
	exit
fi
if [ $CURRENT_MINOR_VERSION -gt $LATEST_MINOR_VERSION ]; then
	echo "Error! Current $CURRENT_VERSION is greater than latest $LATEST_VERSION. Corrupted installation?"
	exit
fi
if [ $CURRENT_MINOR_VERSION -eq $LATEST_MINOR_VERSION ]; then
	echo "Already at latest version, nothing needs to be updated"
	exit
fi


# Start update
echo Will update to version $LATEST_VERSION
INST_PATH=/opt
CHIP_PATH=${INST_PATH}/chipster
TOOLS_PATH=${CHIP_PATH}/tools
TMPDIR_PATH=/tmp/chipster-install-temp

# create temp dir
rm -rf ${TMPDIR_PATH}/
mkdir ${TMPDIR_PATH}/

#
# VERSION SPECIFIC ENTRIES START HERE
# (ADD NEW ENTRIES TO THE END)
#

# 2.0.3
if [ $CURRENT_MINOR_VERSION -lt 3 ] ; then
	echo "Installing mm10 bowtie indexes"
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_mm10.tar.gz  | tar -xz -C ${TOOLS_PATH}/bowtie/
fi

# 2.0.4 
if [ $CURRENT_MINOR_VERSION -lt 4 ] ; then

    echo "Updating prinseq"
		cd ${TMPDIR_PATH}/
    curl -sL http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.19.3.tar.gz/download | tar -xz
    chmod a+x prinseq-lite-0.19.3/prinseq-lite.pl
    chmod a+x prinseq-lite-0.19.3/prinseq-graphs.pl
    mv prinseq-lite-0.19.3 ${TOOLS_PATH}/
    rm ${TOOLS_PATH}/prinseq
    ln -s prinseq-lite-0.19.3 ${TOOLS_PATH}/prinseq
		rm -rf ${TOOLS_PATH}/prinseq-lite-0.17.3

    echo "Adding nz131a520662fcdf"
    cd ${TMPDIR_PATH}/
    wget http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/nz131a520662fcdf.tar.gz
    ${TOOLS_PATH}/R/bin/R CMD INSTALL nz131a520662fcdf.tar.gz
    rm nz131a520662fcdf.tar.gz
    
    echo "Installing vcftools"
		cd ${TMPDIR_PATH}/
		curl -sL http://sourceforge.net/projects/vcftools/files/vcftools_0.1.9.tar.gz/download| tar -xz
		cd vcftools_0.1.9/
		make
		cd ../
		mv vcftools_0.1.9/ ${TOOLS_PATH}/
		ln -s vcftools_0.1.9 ${TOOLS_PATH}/vcftools
    
fi

#
# VERSION SPECIFIC ENTRIES END HERE
#

# Update Chipster itself (incl. tool scripts), unless already at latest
if [ $CURRENT_MINOR_VERSION -lt $LATEST_MINOR_VERSION ] ; then

	echo "Updating Chipster installation to $LATEST_VERSION"

	# Get install package (override, if exists)
	wget -Nq http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/versions/$LATEST_VERSION/chipster-$LATEST_VERSION.tar.gz .

	# Remove old libs to avoid conflicts when lib names change
	rm -rf shared
	rm -rf webstart/web-root/lib

	# Unpack libs
	echo "Updating libs: shared/libs (local changes will be overridden)"
	tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/shared
	echo "Updating libs: webstart/web-root/lib (local changes will be overridden)"
	tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/webstart/web-root/lib

	# Unpack and possibly override tool scripts
	echo "Updating tool scripts: comp/modules (conflicting local changes will be overridden)"
	tar -C .. --overwrite -xzf chipster-$LATEST_VERSION.tar.gz chipster/comp/modules

	# Clean up
	rm chipster-$LATEST_VERSION.tar.gz
fi

# remove temp dir
rm -rf ${TMPDIR_PATH}/

# We are done
echo "Update completed successfully"
