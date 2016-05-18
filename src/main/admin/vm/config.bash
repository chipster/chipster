#!/usr/bin/env bash

#This is config file for chipster installer. 
#It will be loaded at the start of a install-chipster.sh -script

#Logfile
export logfile="../install_logs.txt"

#Module directories
export moduledir="module-list.bash"

#Parallel processing
export parallel="0"

#Parallel jobs
export jobs="4"

#Force install, overrides all existing modules
export force="1"

## Initialize:

# Paths
export CHIP_PATH=/opt/chipster
export TOOLS_PATH=${CHIP_PATH}/tools
# check if TMPDIR is set and set it if necessary
if [ -z ${TMPDIR+x} ]; then
  TMPDIR=/tmp
fi
export TMPDIR
export TMPDIR_PATH=$TMPDIR/install

export CHIPSTER_UID=1001
export CHIPSTER_GID=1001
export UBUNTU_UID=1000
export UBUNTU_GID=1000

# EMBOSS_OPTIONS is used in multiple scripts
export EMBOSS_VERSION=6.5.7
export EMBOSS_PATH=${TOOLS_PATH}/EMBOSS-${EMBOSS_VERSION}
export EMBOSS_OPTIONS="--prefix=${EMBOSS_PATH}"

# nic.funet.fi service endpoint
export NIC_MIRROR=bio.nic.funet.fi
#NIC_MIRROR=www.nic.funet.fi

# tmpdir for parallel
export paralleldir="/mnt/parallel/"
if [ ! -d $paralleldir ]; then
	mkdir $paralleldir
fi



