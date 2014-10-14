#!/usr/bin/env bash

#
# This script will install one module, given in parameter
#
# Notice! This script needs super-user rights!!
# e.g. sudo bash install-module.sh start/chipster_main.bash 2>&1
#

source config.bash
source installation_files/functions.bash

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

# Install module and modules that depend on it
if [ -f $1 ]; then
	install_module $1 2>&1 | tee $logfile
	install_dep $1 2>&1 | tee $logfile
else
echo 'Error: Module ' $1 ' not found.'
fi
