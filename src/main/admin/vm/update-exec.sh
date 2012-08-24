#!/bin/bash

# This script updates to latest minor version of the same major version (e.g. 2.0.1 -> 2.0.3)
LATEST_VERSION=2.0.3

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

# 2.0.3: Updated mouse genome
if [ $CURRENT_MINOR_VERSION -lt 3 ] ; then
	echo "Installing mm10 bowtie indexes"
	curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/bowtie_indexes/bowtie_index_mm10.tar.gz  | tar -xz -C ${TOOLS_PATH}/bowtie/
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

	# Unpack libs
	echo "Updating libs: shared/libs (local changes will be overridden)"
	tar -C .. -xzf chipster-2.0.3.tar.gz chipster/shared

	# Unpack and possibly override tool scripts
	echo "Updating tool scripts: comp/modules (conflicting local changes will be overridden)"
	tar -C .. --overwrite -xzf chipster-2.0.3.tar.gz chipster/comp/modules
	
	# Update client libs
	# TODODODODODO

	# Clean up
	rm chipster-$LATEST_VERSION.tar.gz
fi

# remove temp dir
rm -rf ${TMPDIR_PATH}/

# We are done
echo "Update completed successfully"
