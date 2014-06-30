#!/usr/bin/env bash

#
# This script will install Chipster 3, w/ dependencies
#
# Notice! This script needs super-user rights!!
# e.g. sudo bash install-chipster.sh 2>&1
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


## If parallel is ON, install GNU parallel
if [ $parallel == "1" ]; then
	bash installation_files/install_parallel.bash 2>&1 | tee $logfile
fi


echo "Starting to generate installation files"
source installation_files/generate_installation_file.bash start/ inst_files/start.bash
source installation_files/generate_installation_file.bash external_genomes/ inst_files/gnomes.bash
source installation_files/generate_installation_file.bash external_indexes/ inst_files/indx.bash
source installation_files/generate_installation_file.bash external_tools/ inst_files/tools.bash
source installation_files/generate_installation_file.bash R/ inst_files/R.bash
source installation_files/generate_installation_file.bash finish/ inst_files/finish.bash

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


# Install start modules
install_file inst_files/start.bash 2>&1 | tee $logfile


#Install R with libraries
install_file inst_files/R.bash 2>&1 | tee $logfile


#install tools
install_file inst_files/tools.bash 2>&1 | tee $logfile


#Make symbolic links
bash ../installation_files/sym_link_to_admin_scripts.bash 2>&1 | tee $logfile


#External genomes
install_file inst_files/gnomes.bash 2>&1 | tee $logfile


#External indeces
install_file inst_files/indxs.bash 2>&1 | tee $logfile


#Wait for R libs -installation
wait

#Finish installation
install_file inst_files/finish.bash 2>&1 | tee $logfile

# Return to installation directory
cd ..

echo "Installation completed"

