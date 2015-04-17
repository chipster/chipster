#!/bin/bash

# first set workin directory to the directory of the script
cd "$(dirname "$0")"

# paths in configuration file assume parent working directory
cd ..

# if there are parameters
if (($#))
then
	# set config path
	CONF="--config=conf/chipster-config.xml"
fi

# command for running java
COMMAND="sudo -u chipster java -cp ../shared/lib/*: fi.csc.microarray.messaging.admin.CompAdmin"

# run and hide help texts of config parameter, because it's already set by this script
$COMMAND $CONF $@ | sed "s/--config=CHIPSTER-CONFIG.XML//g" | sed "s/Chipster manager config file//g"
