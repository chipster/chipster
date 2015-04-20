#!/bin/bash

# This script updates to latest version. Updates between minor versions should be smooth and
# automatic, where as updates between major versions can require some manual steps afterwards
# if some specific local customizations were in place.
# This update mechanism has been available since 2.0.2.

# Latest version, matching tar-packages must be available 
##
LATEST_VERSION=3.4.0
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

# Check current version
echo Detected version $CURRENT_VERSION
compare_to_current "$LATEST_VERSION"
if [ $CURRENT_COMPARED -gt 0 ] ; then 
  echo "Update error: current version $CURRENT_VERSION is newer than latest $LATEST_VERSION"
  exit 1
fi
if [ $CURRENT_COMPARED -eq 0 ] ; then 
  echo "Already at the latest version, nothing needs to be updated"
  exit
fi
echo "Will update to version $LATEST_VERSION"
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
chmod a+rwx ${BACKUPDIR_PATH} # both ubuntu and chipster need to be able to write here

#######################################
# VERSION SPECIFIC ENTRIES START HERE #
# (ADD NEW ENTRIES TO THE END)        #
#######################################


# 3.1.0
compare_to_current_and_latest "3.1.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then
  echo ""
fi

# 3.4.0
compare_to_current_and_latest "3.4.0"
if [ $CURRENT_COMPARED -lt 0 ] && [ ! $LATEST_COMPARED -lt 0 ] ; then

    echo ""
    echo "Download the latest virtual machine"
    echo ""
    echo "Because of the following changes, an automatic update isn't available for this version: "
    echo "- update genomes in tools package to latest Ensembl release"
    echo "- install jobmanager to support long running jobs"
    echo ""
    echo "Please download the latest virtual machine from http://chipster.github.io/chipster/"
    echo ""
    exit 0
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
  sudo rm -f chipster-$LATEST_VERSION.tar.gz
  sudo -u chipster wget http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/versions/$LATEST_VERSION/chipster-$LATEST_VERSION.tar.gz

  # Unpack libs
  echo "** Updating Chipster libs: shared/libs"
  sudo -u chipster mv shared ${BACKUPDIR_PATH}/
  sudo -u chipster tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/shared

  # Unpack webstat web-root including client jar
  echo "** Updating Chipster web: webstart/web-root"
  sudo -u chipster cp webstart/web-root/chipster.jnlp ${BACKUPDIR_PATH}/
  sudo -u chipster cp webstart/web-root/chipster-config.xml ${BACKUPDIR_PATH}/
  sudo -u chipster mv webstart/web-root ${BACKUPDIR_PATH}/
  sudo -u chipster tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/webstart/web-root
  sudo -u chipster cp ${BACKUPDIR_PATH}/chipster.jnlp webstart/web-root/ 
  sudo -u chipster cp ${BACKUPDIR_PATH}/chipster-config.xml webstart/web-root/ 

  # Copy away tool scripts in case there were important local changes
  sudo -u chipster cp -r comp/modules ${BACKUPDIR_PATH}/

  # Unpack tool scripts
  echo "** Updating Chipster tool scripts: comp/modules"
  sudo -u chipster tar -C .. --overwrite -xzf chipster-$LATEST_VERSION.tar.gz chipster/comp/modules

  # Update runtimes.xml
  echo "** Updating Chipster runtimes: comp/conf/runtimes.xml"
  sudo -u chipster cp -r comp/conf/runtimes.xml ${BACKUPDIR_PATH}/
  sudo -u chipster tar -C .. --overwrite -xzf chipster-$LATEST_VERSION.tar.gz chipster/comp/conf/runtimes.xml

  # Update webapps
  sudo -u chipster rm -rf webstart/webapps
  sudo -u chipster tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/webstart/webapps/tool-editor.war
  sudo -u chipster rm -rf manager/webapps
  sudo -u chipster tar -C .. -xzf chipster-$LATEST_VERSION.tar.gz chipster/manager/webapps/admin-web.war

  # Clean up
  sudo -u chipster rm chipster-$LATEST_VERSION.tar.gz
fi

# Remove temp dir
rm -rf ${TMPDIR_PATH}

# Check backup dir
SIZE=`du -hs ${BACKUPDIR_PATH} | cut -f1`
echo ""
echo "Total of $SIZE old data has been backed up to ${BACKUPDIR_PATH}"
echo "It is recommended to inspect the directory and then to remove it"
   
# We are done
echo "Update completed successfully"
echo "Remember to start the Chipster service: 'service chipster start'"
echo $END_MESSAGE
