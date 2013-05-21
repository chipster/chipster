#!/bin/bash

# This script updates to latest version. Updates between minor versions should be smooth and
# automatic, where as updates between major versions can require some manual steps afterwards
# if some specific local customisations were in place.
# This update mechanism has been available since 2.0.2.

# Update to version
LATEST_VERSION=2.5.2

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
# Check version from somewhere
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


# 2.5.2
compare_to_current "2.5.2"
if [ $CURRENT_COMPARED -lt 0 ] ; then 

  echo "** Updating HTSeq"
  pip install HTSeq==0.5.4p3
  wget -O /usr/local/bin/htseq-count_chr http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/htseq/htseq-count_chr 
  chmod 755 /usr/local/bin/htseq-count_chr
  wget -O /usr/local/lib/python2.7/dist-packages/HTSeq/scripts/count_chr.py http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/htseq/count_chr.py

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
