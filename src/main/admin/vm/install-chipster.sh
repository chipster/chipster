#!/usr/bin/env bash

#
# This script will install Chipster 3, w/ dependencies
#
# Notice! This script needs super-user rights!!
# e.g. sudo bash install-chipster.sh
#



# Load config-file
source config.bash


# Load installation functions
source installation_files/functions.bash


## Create tmpdir
rm -rf ${TMPDIR_PATH}
mkdir -p ${TMPDIR_PATH}/


## Create installation file directories
rm -rf modules/inst_files/
mkdir modules/inst_files/

## Create flag-directories
create_flag_dirs

## Load module directories to array
read_dirs foldiers.bash


## If parallel is ON, install GNU parallel
if [ "$parallel" == "1" ]; then
	bash installation_files/install_parallel.bash 2>&1 | tee $logfile
fi


echo "Starting to generate installation files"

for i in "${foldiers[@]}"; do
	source installation_files/generate_installation_file.bash $i/ inst_files/$i.bash
done

echo "Installation files rdy, proceed to installation"


# Set execution trace
set -x

# Set exit on error
set -e

# Set fail on pipe
set -o pipefail

# Set exit on unset variable
set -u

# Move to modules working directory
cd modules/


# Install modules

for i in "${foldiers[@]}"; do
	
	# Wait for R
	if [ "$i" == "finish" ]; then
		wait
	fi
	install_file inst_files/$i.bash 2>&1 | tee $logfile
done


# Return to installation directory
cd ..

echo "Installation completed"

