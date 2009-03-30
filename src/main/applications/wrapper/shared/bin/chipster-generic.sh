#!/bin/bash

CHIPSTER_HOME="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/../.."
COMPONENT=$1
COMMAND=$2

# platform, should be linux-x86-32 or linux-x86-64, autodetected if not set
PLATFORM=""


## detect if platform not set
if [ ! $PLATFORM ]; then
    ARCH=`uname -m`
    if [ "$ARCH" = "i686" ]; then
	PLATFORM="linux-x86-32"
    elif [ "$ARCH" = "x86_64" ]; then
	PLATFORM="linux-x86-64"
    else
	echo "Could not detect hardware architecture, please set platform manually."
	exit -1;
    fi
fi


# run command
$CHIPSTER_HOME/$COMPONENT/bin/$PLATFORM/chipster-$COMPONENT $COMMAND

