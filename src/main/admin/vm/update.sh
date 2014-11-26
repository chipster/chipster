#!/bin/bash

# Bootstrap script for updating Chipster installation.
# This script only downloads the actual update script, allowing the update mechanism to be updated also.

echo ""
echo "Some parts of the update may need root privileges. These parts are run using sudo."
echo "Testing permission to use sudo..."
if [ "$(sudo whoami)" != 'root' ]; then echo 'You need sudo rights to run the update script, aborting.'; 
exit 1; fi
echo "Sudo ok"
echo ""

# Update file web location and name
UPDATE_URL_PREFIX=http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/virtual_machines/updates
MAIN_UPDATE_FILE=update-exec3.sh

# Figure out Chipster major version and corresponding URL
CHIPSTER_VERSION=`ls -1 shared/lib | grep ^chipster-[0-9\\.]*.jar | gawk 'match($0, "chipster-([0-9\\\\.]*).jar", g) {print g[1]}'`
CHIPSTER_MAJOR_VERSION=`echo $CHIPSTER_VERSION | gawk 'match($0, "([0-9]*.[0-9]*).[0-9]*", g) {print g[1]}'`
MAIN_UPDATE_FILE_URL=$UPDATE_URL_PREFIX/$CHIPSTER_MAJOR_VERSION/$MAIN_UPDATE_FILE

# Remove old update file, if exists
if [ -e $MAIN_UPDATE_FILE ]; then
  sudo rm $MAIN_UPDATE_FILE
fi

# Download latest update file
sudo wget -q $MAIN_UPDATE_FILE_URL

# Check if download was ok
if [ -e $MAIN_UPDATE_FILE ]; then

	# Run update
	sudo chmod a+x $MAIN_UPDATE_FILE
	sudo -u ubuntu ./$MAIN_UPDATE_FILE
	
	# Clean up
	sudo rm $MAIN_UPDATE_FILE
	
	exit 0

else

	echo Failed to download latest update file from $MAIN_UPDATE_FILE_URL
	exit 1
fi